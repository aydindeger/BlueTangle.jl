# Independent binary error propagation for the actual Shor extraction schedule.
# This deliberately does not call propagate_pauli, pauli_syndrome, or the
# production decoder. Individual cat readouts are random, but their parity and
# the verifier result are deterministic for a Pauli error on a code state.
@testset "Circuit-level single-fault audit of Shor recovery" begin
    code = steane_code()
    circuits = [shor_syndrome_circuit(code, i) for i in 1:code.m]
    n = code.n + 5
    paulis = ("I", "X", "Y", "Z")
    H = code.generator[1:3, 1:code.n]

    function add_paulis!(x, z, qubits, names)
        for (q, name) in zip(qubits, names)
            name in ("X", "Y") && (x[q] ⊻= 1)
            name in ("Z", "Y") && (z[q] ⊻= 1)
        end
    end
    wires(op) = op.q == 1 ? (op.qubit,) : (op.qubit, op.target_qubit)
    function conjugate!(x, z, op)
        q, t = op.qubit, op.target_qubit
        if op.name == "H"
            x[q], z[q] = z[q], x[q]
        elseif op.name in ("CX", "CNOT")
            x[t] ⊻= x[q]
            z[q] ⊻= z[t]
        elseif op.name == "CZ"
            z[q] ⊻= x[t]
            z[t] ⊻= x[q]
        elseif op.name in paulis
            # A Pauli gate leaves an error frame unchanged up to phase.
            nothing
        else
            error("audit needs an independent rule for $(op.name)")
        end
    end
    function syndrome(x, z)
        vcat(mod.(H * z[1:7], 2), mod.(H * x[1:7], 2))
    end
    function recovery(bits)
        sx, sz = bits[1:3], bits[4:6]
        qx = all(iszero, sz) ? 0 : findfirst(q -> H[:, q] == sz, 1:7)
        qz = all(iszero, sx) ? 0 : findfirst(q -> H[:, q] == sx, 1:7)
        [Op(q == qx ? (q == qz ? "Y" : "X") : "Z", q)
         for q in 1:7 if q == qx || q == qz]
    end
    function corrected_logical(x, z)
        x, z = copy(x), copy(z)
        for op in recovery(syndrome(x, z))
            add_paulis!(x, z, (op.qubit,), (op.name,))
        end
        # In this encoder basis, logical X and Z both have support on 5,6,7.
        return (x=sum(x[5:7]) % 2, z=sum(z[5:7]) % 2)
    end
    function frame_run(; fault_location=0, fault_names=(), initial=Op[])
        x, z = zeros(Int, n), zeros(Int, n)
        for op in initial
            add_paulis!(x, z, (op.qubit,), (op.name,))
        end
        schedule = NamedTuple[]
        history = Vector{Int}[]
        rejections = 0
        function execute(event, round, check, attempt)
            op, kind = event.op, event.kind
            push!(schedule, (; event..., round, check, attempt))
            location = length(schedule)
            if kind == :prep
                x[op.qubit] = z[op.qubit] = 0
            elseif kind == :gate
                conjugate!(x, z, op)
            elseif kind != :measure
                error("unknown location kind: $kind")
            end
            # Measurements carry a fault immediately before their outcome;
            # preparation and gate faults occur immediately after the operation.
            location == fault_location && add_paulis!(x, z, wires(op), fault_names)
            return kind == :measure ? x[op.qubit] : 0
        end
        for round in 1:6
            bits = Int[]
            for check in 1:code.m
                events = circuits[check]
                verify = findfirst(e -> e.kind == :measure && e.stage == :cat_verify, events)
                @assert verify !== nothing
                attempt = 0
                while true
                    attempt += 1
                    verifier = 0
                    for event in events[1:verify]
                        verifier = execute(event, round, check, attempt)
                    end
                    verifier == 0 && break
                    rejections += 1
                end
                bit = 0
                for event in events[verify+1:end]
                    outcome = execute(event, round, check, attempt)
                    event.kind == :measure && (bit ⊻= outcome)
                end
                push!(bits, bit)
            end
            push!(history, bits)
            if length(history) >= 2 && bits == history[end-1]
                for op in recovery(bits)
                    # The production recovery is a physical Pauli gate. First
                    # update the reference error, then add a possible gate fault.
                    add_paulis!(x, z, (op.qubit,), (op.name,))
                    execute((kind=:gate, stage=:recovery, op=op), round, 0, 0)
                end
                return (; x=x[1:7], z=z[1:7], schedule, history, rejections)
            end
        end
        error("a single-fault trajectory failed to reach agreement")
    end

    baseline = frame_run()
    @test length(baseline.schedule) == 288
    @test length(baseline.history) == 2
    fault_count = rejected_count = residual_count = max_rounds = 0
    for (location, event) in enumerate(baseline.schedule)
        for names in Iterators.product(ntuple(_ -> paulis, length(wires(event.op)))...)
            all(==("I"), names) && continue
            result = frame_run(; fault_location=location, fault_names=names)
            fault_count += 1
            rejected_count += result.rejections > 0
            residual_count += any(!iszero, syndrome(result.x, result.z))
            max_rounds = max(max_rounds, length(result.history))
            @test corrected_logical(result.x, result.z) == (x=0, z=0)
        end
    end
    @test fault_count == 2160
    @test rejected_count == 552
    @test residual_count == 102
    @test max_rounds == 4

    # Check the real statevector implementation against representative frame
    # predictions, including rejected preparation, readout, coupling, and a
    # late fault leaving a correctable residual. Only syndrome parities are
    # compared: individual cat outcomes are intrinsically random.
    logical_input = normalize(ComplexF64[1, 0.2+0.7im])
    encoded = code.encode(logical_input)
    for (location, names) in ((8, ("X", "I")), (12, ("X",)),
                               (13, ("X", "Y")), (23, ("X",)),
                               (280, ("I", "Y")))
        predicted = frame_run(; fault_location=location, fault_names=names)
        fault = loc -> loc.index == location ?
            [Op(name, q; noisy=false) for (q, name) in zip(wires(loc.op), names) if name != "I"] : Op[]
        actual = shor_qec_round(code, encoded; fault=fault)
        @test actual.history == predicted.history
        @test sum(sum(attempts .- 1) for attempts in actual.attempts) == predicted.rejections
        @test length(actual.trace) == length(predicted.schedule)
        diagnosed = qec_round(code, actual.state; decoder=steane_decoder(code)).state
        @test abs2(dot(logical_input, code.decode(diagnosed))) ≈ 1
    end

    # The complementary distance-three EC condition: one incoming data error
    # and no internal fault is removed completely, including its logical action.
    for q in 1:7, name in ("X", "Y", "Z")
        result = frame_run(; initial=[Op(name, q)])
        @test all(iszero, syndrome(result.x, result.z))
        @test corrected_logical(result.x, result.z) == (x=0, z=0)
    end

    # Product |0> and |+> preparation projects onto one of eight initial check
    # sectors. Audit all sectors rather than silently assuming the zero record.
    # A logical Z is harmless on |0>, and a logical X is harmless on |+>.
    for basis in (:zero, :plus)
        actual_preparation = steane_prepare(basis)
        target = basis == :zero ? zero_state(1) : plus_state(1)
        @test abs2(dot(target, code.decode(actual_preparation.state))) ≈ 1
        preparation_events = filter(e -> e.stage in (:data_prep, :data_basis), actual_preparation.trace)
        @test length(preparation_events) == (basis == :zero ? 7 : 14)
        for sector in 0:7
            sector_bits = [Int((sector >> (3-i)) & 1) for i in 1:3]
            bits = basis == :zero ? vcat(sector_bits, zeros(Int, 3)) :
                                   vcat(zeros(Int, 3), sector_bits)
            initial = recovery(bits)
            prepared = frame_run(; initial)
            for (location, event) in enumerate(prepared.schedule)
                for names in Iterators.product(ntuple(_ -> paulis, length(wires(event.op)))...)
                    all(==("I"), names) && continue
                    result = frame_run(; initial, fault_location=location, fault_names=names)
                    logical = corrected_logical(result.x, result.z)
                    @test (basis == :zero ? logical.x : logical.z) == 0
                    @test length(result.history) <= 4
                end
            end
            # Faults during physical product preparation and its optional H
            # layer must also be propagated, rather than assuming ideal input.
            for (location, event) in enumerate(preparation_events), name in ("X", "Y", "Z")
                x, z = zeros(Int, 7), zeros(Int, 7)
                add_paulis!(x, z, wires(event.op), (name,))
                for later in preparation_events[location+1:end]
                    if later.kind == :prep
                        x[later.op.qubit] = z[later.op.qubit] = 0
                    else
                        conjugate!(x, z, later.op)
                    end
                end
                errors = [Op(x[q] == 1 ? (z[q] == 1 ? "Y" : "X") : "Z", q)
                          for q in 1:7 if x[q] != 0 || z[q] != 0]
                result = frame_run(; initial=vcat(initial, errors))
                logical = corrected_logical(result.x, result.z)
                @test (basis == :zero ? logical.x : logical.z) == 0
                @test all(iszero, syndrome(result.x, result.z))
            end
        end
    end
end
