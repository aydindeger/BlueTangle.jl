@testset "Independent algebra audit of the encoded factory" begin
    # The production helper gives only the CNOT schedule. This test applies it
    # with independent bit permutations and compares it with the diagonal
    # operator obtained directly from the fifteen Reed–Muller column labels.
    steps = BlueTangle._steane_distillation15_steps()
    @test length(steps) == 16
    @test sort([s.column for s in steps[1:15]]) == collect(1:15)
    @test last(steps).column == 0
    @test sum(length(s.controls) for s in steps) == 16

    function independent_cnot(psi, control)
        result = similar(psi)
        for label in 0:31
            # The four check wires are MSB-first; output wire 5 is the low bit.
            destination = label ⊻ (((label >> (5-control)) & 1) == 1 ? 1 : 0)
            result[destination+1] = psi[label+1]
        end
        return result
    end
    function diagonal_output(psi, diagonal)
        [diagonal[(label & 1)+1] * psi[label+1] for label in 0:31]
    end
    function run_schedule(psi, theta; bad_column=0, resource=ComplexF64[1, 0], branch=0)
        phase = cis(theta)
        for step in steps
            for control in step.controls
                psi = independent_cnot(psi, control)
            end
            step.column == 0 && continue
            # Resource readout 0 gives diag(a,b); readout 1, after the S/SD
            # feedback, gives diag(b, exp(2i*theta)*a). No normalization is
            # applied: retained amplitudes carry their acceptance probabilities.
            diagonal = step.column == bad_column ?
                (branch == 0 ? resource : [resource[2], phase^2 * resource[1]]) : [1, phase]
            psi = diagonal_output(psi, diagonal)
        end
        return diagonal_output(psi, [1, phase^2])
    end
    project_checks(psi) = ComplexF64[sum(psi[2label+bit+1] for label in 0:15) / 4 for bit in 0:1]

    plus5 = fill(ComplexF64(inv(sqrt(32))), 32)
    for theta in (pi/4, -pi/4)
        independent_diagonal = ComplexF64[]
        for label in 0:31
            angle_count = 0
            for column in 1:15
                parity = label & 1
                for q in 1:4
                    parity ⊻= ((label >> (5-q)) & 1) & ((column >> (4-q)) & 1)
                end
                angle_count += parity
            end
            # Final S (or SD) supplies the last two eighth-turns on output.
            push!(independent_diagonal, cis(theta * (angle_count + 2(label & 1))))
        end
        basis = Matrix{ComplexF64}(I, 32, 32)
        implemented = hcat([run_schedule(basis[:, j], theta) for j in 1:32]...)
        @test isapprox(implemented, Diagonal(independent_diagonal); atol=1e-12)
        A = ComplexF64[1, cis(theta)] / sqrt(2)
        @test isapprox(project_checks(run_schedule(plus5, theta)), A; atol=1e-12)

        # Computational resource vectors span every pure/mixed one-qubit
        # resource defect. In either feed-forward branch, the accepted map has
        # range span(A). Their summed effect is the A-projector: an orthogonal
        # resource defect is rejected, and any accepted output is the right A.
        for column in 1:15
            maps = [hcat([project_checks(run_schedule(plus5, theta;
                        bad_column=column, resource=resource, branch=branch))
                         for resource in (ComplexF64[1, 0], ComplexF64[0, 1])]...)
                    for branch in 0:1]
            for K in maps
                @test norm((I - A*A') * K) < 1e-12
            end
            @test isapprox(sum(K' * K for K in maps), A*A'; atol=1e-12)
        end
    end
end

@testset "Raw encoder single faults become logical-resource defects" begin
    code = steane_code()
    q = only(code.logical_indices)
    schedule = vcat([Op("I", wire; noisy=false) for wire in 1:7],
                    [Op("H", q), Op("T", q)], code.ops_encoding)
    before = [zero_state(7)]
    for op in schedule
        push!(before, apply(last(before), op))
    end
    A = ComplexF64[1, cis(pi/4)] / sqrt(2)
    @test abs2(dot(A, code.decode(last(before)))) ≈ 1
    fault_count = logical_defects = 0
    for (location, op) in enumerate(schedule)
        wires = op.q == 1 ? (op.qubit,) : (op.qubit, op.target_qubit)
        for names in Iterators.product(ntuple(_ -> ("I", "X", "Y", "Z"), length(wires))...)
            all(==("I"), names) && continue
            faults = [Op(name, wire; noisy=false) for (wire, name) in zip(wires, names) if name != "I"]
            faulty = apply(before[location+1], faults)
            for later in schedule[location+1:end]
                faulty = apply(faulty, later)
            end
            # Use the real Shor circuit with no additional faulty locations:
            # this is the single-fault condition for an error in raw encoding.
            corrected = shor_qec_round(code, faulty).state
            logical = code.decode(corrected)
            @test isapprox(norm(logical), 1; atol=1e-10)
            @test isapprox(code.encode(logical), corrected; atol=1e-10)
            fault_count += 1
            logical_defects += abs2(dot(A, logical)) < 1 - 1e-8
        end
    end
    @test fault_count > 100
    @test logical_defects > 0 # raw preparation is deliberately not an FT factory
end

@testset "Factory trace preserves component and aborted work" begin
    trace = NamedTuple[]
    tags = (factory_attempt=1, component=:raw, block=0, resource=3)
    observed = NamedTuple[]
    child = BlueTangle._factory_call(trace, tags, event -> (push!(observed, event); Op[])) do callback
        raw_magic_state(; phase_error=1.0, fault=callback)
    end
    @test length(trace) == length(observed) == length(child.trace) == 3
    @test getproperty.(trace, :index) == [1, 2, 3]
    @test all(event -> event.resource == 3 && event.component == :raw, trace)
    @test trace[end].phase_error_applied
    @test trace[end].outcome === nothing

    # Record a second component before a heralded interruption. Previously
    # completed components keep their records and global indices stay unique.
    @test_throws QECProtocolFailure BlueTangle._factory_call(trace, tags, nothing) do callback
        callback((index=1, kind=:measure, stage=:test, round=0, check=0,
                  attempt=1, op=Op("MZ", 1)))
        throw(QECProtocolFailure("deliberate test interruption"))
    end
    @test getproperty.(trace, :index) == [1, 2, 3, 4]
    @test trace[end].local_index == 1
    @test trace[end].outcome === nothing
end
