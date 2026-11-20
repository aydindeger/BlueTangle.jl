@testset "Backend integration and MPS ordering" begin
    rng = MersenneTwister(92)
    sites = N_MPS(4)
    vector = normalize(randn(rng, ComplexF64, 16))
    psi = to_MPS(vector, sites; cutoff=0.0)
    @test to_state(psi) ≈ vector
    @test inner(psi, im*vector) ≈ im
    @test inner(im*vector, psi) ≈ -im
    @test fidelity(psi,psi) ≈ 1
    @test fidelity(psi, vector) ≈ 1
    @test fidelity(vector, psi) ≈ 1
    @test correlation(psi,"X,Y",[1,4]) ≈ correlation(vector,"X,Y",[1,4])
    @test expect(psi,Op("Z",2)) ≈ expect(vector,Op("Z",2))
    before = to_state(psi)
    @test length(sample_bit(psi, 5)) == 5
    @test to_state(psi) ≈ before
    @test entanglement_entropy(zero_state(N_MPS(1))) == 0
    @test sample_exact(psi)[2] ≈ abs2.(vector)
    @test measure(psi).number_of_experiment == -1
    @test measure(psi).sample ≈ abs2.(vector)

    large = one_state(N_MPS(72))
    measurement = measure(large, 3)
    @test measurement.bitstr == [big(2)^72-1]
    @test measurement.expect == -ones(72)
    @test isempty(measurement.mid_measurements)
    @test_throws ArgumentError sample_exact(large)
    @test_throws ArgumentError measure(large)
    @test int2bin(bin2int(ones(Int,72)),72) == ones(Int,72)

    for backend in (:statevector, :tensor, :density_matrix, :camps)
        @test measure([Op("I",1)], 20; backend=backend,
                      options=Options(readout_noise=Noise1("bit_flip",1.0))).bitstr == [1]
        @test measure([Op("X",1),Op("MZ",1),Op("RES",1),Op("MZ",1)],2;backend=backend).mid_measurements == [[1,0],[1,0]]
        reset = [Op("H",1),Op("CX",1,2),Op("RES",1)]
        result = measure(reset, 1000; backend=backend)
        @test result.expect[1] ≈ 1
        @test abs(result.expect[2]) < 0.15
        @test isempty(result.mid_measurements)
    end
    # A Kraus channel must be sampled afresh on each trajectory, even when no
    # NoiseModel option or measurement gate is present in the circuit.
    noisy = [Op("I",1),OpQC("bit_flip",0.5,1)]
    @test abs(measure(noisy,1000).expect[1]) < 0.15
    @test abs(measure(noisy,1000;backend=:density_matrix).expect[1]) < 0.15

    mixed = ComplexF64[0.25 0; 0 0.75]
    @test measure(mixed).sample ≈ [0.25,0.75]
    @test measure(mixed).expect ≈ [-0.5]
    @test (Op("X",1) * mixed) ≈ gate.X * mixed * gate.X

    options = Options(readout_noise=Noise1("bit_flip",0.1),measurement_mitigate=true)
    mitigated = measure([Op("I",1)], 1; options=options)
    @test length(mitigated.sample) == 2
    @test sum(mitigated.sample) ≈ 1
    @test minimum(mitigated.sample) < 0
    @test abs(only(mitigated.expect)) ≈ 1.25
    # A finite-shot histogram need not lie in the image of physical input
    # probabilities. Inversion must retain signed weights to undo the response.
    _, signed = BlueTangle._mitigate_probabilities([0],[1.0],1,options.readout_noise)
    response = [0.9 0.1; 0.1 0.9]
    @test signed ≈ [1.125,-0.125]
    @test response * signed ≈ [1.0,0.0] atol=1e-14
    _, signed_pair = BlueTangle._mitigate_probabilities([0],[1.0],2,options.readout_noise)
    @test kron(response,response) * signed_pair ≈ [1.0,0.0,0.0,0.0] atol=1e-14
    bits, corrected = BlueTangle._mitigate_probabilities([0,1,2,3],[0.81,0.09,0.09,0.01],2,options.readout_noise)
    @test bits == 0:3
    @test corrected ≈ [1,0,0,0] atol=1e-14

    for op in (Op("custom",gate.CX,1,3), Op("custom",gate.CY,3,1), Op("SWAP",1,3;control=2))
        @test sa.sparse([op]) ≈ op.expand(3)
    end
    conditional = ifOp("MZ",1,Op("custom",gate.CX,1,3),Op("I",1))
    ops = [Op("H",3),conditional]
    @test apply(zero_state(3), ops) ≈ kron(zero_state(2), plus_state(1))
    @test isempty([name for name in names(BlueTangle) if !isdefined(BlueTangle,name)])
    @test :run ∉ names(BlueTangle)
    if isempty(methods(phasepoly))
        @test_throws ArgumentError apply([Op("H", 1)]; backend=:phasepoly)
        @test_throws ArgumentError measure([Op("H", 1)], 10; backend=:phasepoly)
        @test_throws ArgumentError sample([Op("H", 1)], 10; backend=:phasepoly)
        @test_throws ArgumentError shadow([Op("H", 1)], 10; backend=:phasepoly)
    end
end

@testset "Classical readout mitigation requirements" begin
    # Coherent rotations can change measured populations through off-diagonal
    # entries of rho; a classical confusion matrix cannot undo that dependence.
    for model in ("rot_x", "rot_y")
        channel = Noise1(model, pi/3)
        @test_throws ArgumentError BlueTangle._mitigate_probabilities([0, 1], [0.5, 0.5], 1, channel)
        @test_throws ArgumentError BlueTangle._measurement_mitigate_inv_matrix(1, channel)
    end
    options = Options(readout_noise=Noise1("rot_y", pi/3), measurement_mitigate=true)
    @test_throws ArgumentError measure([Op("H", 1)], 1; options, backend=:density_matrix)

    # Completely random flips and complete relaxation erase the input bit.
    for channel in (Noise1("bit_flip", 0.5), Noise1("amplitude_damping", 1.0))
        @test_throws ArgumentError BlueTangle._mitigate_probabilities([0, 1], [0.5, 0.5], 1, channel)
        @test_throws ArgumentError BlueTangle._measurement_mitigate_inv_matrix(1, channel)
    end

    # Opposite coherent rotations have cancelling coherence contributions.
    # Validity must be checked on the whole channel, not each Kraus operator.
    symmetric = custom_noise(1, "symmetric rotations", [gates("RY(0.4)")/sqrt(2),
                                                        gates("RY(-0.4)")/sqrt(2)])
    for channel in (Noise1("bit_flip", 0.1), Noise1("amplitude_damping", 0.3),
                    Noise1("rot_z", 0.4), symmetric)
        state = normalize(ComplexF64[1, 2+im])
        rho = state * state'
        observed = real.(diag(channel.apply(rho, 1)))
        _, corrected = BlueTangle._mitigate_probabilities([0, 1], observed, 1, channel)
        @test corrected ≈ real.(diag(rho))
        @test BlueTangle._measurement_mitigate_inv_matrix(1, channel) * observed ≈ corrected
    end
end

@testset "Direct operation-list API" begin
    for removed in (:compile, :Circuit, :quantum_circuit)
        @test !isdefined(BlueTangle, removed)
    end
    @test isempty(filter(pair -> any(method -> method.name == :apply, pair),
                         Test.detect_ambiguities(BlueTangle)))

    # A nonlocal entangling gate acts directly on its named qubits. The same
    # operation list works with a vector, a density matrix, or an MPS.
    operations = [Op("H", 1), Op("CX", 1, 3)]
    expected = ComplexF64[1, 0, 0, 0, 0, 1, 0, 0] / sqrt(2)
    initial = zero_state(3)
    @test apply(initial, operations) ≈ expected
    @test apply(initial * initial', operations) ≈ expected * expected'
    @test to_rho(expected) ≈ expected * expected'
    @test to_state(apply(zero_state(N_MPS(3)), operations; cutoff=0.0)) ≈ expected
    @test initial == zero_state(3)
    @test get_N(operations) == 3
    @test length(operations) == 2

    # Routing is an explicit transformation; apply and sampling accept it as an option.
    layout = Layout(ones(Int, 1, 3))
    routed = route(operations, layout)
    @test length(routed) > length(operations)
    @test apply(initial, routed) ≈ expected
    @test apply(operations; backend=:statevector, layout=layout) ≈ expected
    @test apply(operations; backend=:density_matrix, layout=layout) ≈ expected * expected'
    @test to_state(apply(operations; backend=:tensor, layout=layout, cutoff=0.0)) ≈ expected
    @test fidelity(apply(operations; backend=:camps, layout=layout), expected) ≈ 1
    @test length(operations) == 2

    # Literal list order matters: layer scheduling could move MZ(2) before
    # MZ(1), reversing the recorded bits for this program.
    measured = [Op("X", 1), Op("MZ", 1), Op("MZ", 2)]
    for backend in (:statevector, :density_matrix, :mps, :camps)
        @test measure(measured, 3; backend=backend).mid_measurements == fill([1, 0], 3)
        @test measure([Op("X", 1)], 2; N=3, backend=backend).bitstr == [4]
        @test measure(Op[], 2; N=2, backend=backend).bitstr == [0]
        _, bits = apply(measured; backend=backend, track_measurements=true)
        @test bits == [1, 0]
    end
    for state in (zero_state(2), zero_state(2; backend=:density_matrix),
                  zero_state(2; backend=:mps), zero_state(2; backend=:camps))
        _, bits = apply(state, measured; track_measurements=true)
        @test bits == [1, 0]
    end
    damping = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    @test apply(zero_state(1), [Op("X", 1)]; noise=damping) ≈ zero_state(1)
    @test apply(ComplexF64[1 0; 0 0], [Op("X", 1)]; noise=damping) ≈ ComplexF64[1 0; 0 0]
    @test to_state(apply(zero_state(N_MPS(1)), [Op("X", 1)]; noise=damping)) ≈ zero_state(1)
    @test apply(to_rho(one_state(1)), [Op("I", 1)];
                noise=NoiseModel(Noise1("amplitude_damping", 0.3), false)) ≈ ComplexF64[0.3 0; 0 0.7]

    # The existing full-state shadow estimator also accepts operations directly.
    snapshot = shadow([Op("I", 1)], 1)
    @test snapshot isa sa.SparseMatrixCSC
    @test eigvals(Hermitian(Matrix(snapshot))) ≈ [-1, 2]
end

@testset "Common operation constructors" begin
    conditional = Op("MZ",1;if0="I",if1="X")
    @test conditional isa ifOp
    @test apply(one_state(1),conditional) ≈ zero_state(1)
    @test Op(Noise1("bit_flip",1.0),1) isa OpQC
    @test apply(zero_state(1),Op(Noise1("bit_flip",1.0),1)) ≈ one_state(1)
    block = Op("prepare",[Op("X",3)])
    @test block isa OpF
    @test get_N([block]) == 3
    @test apply(zero_state(3), [block]) ≈ apply(zero_state(3),Op("X",3))
    @test apply(zero_state(1),Op("full gate",gate.X)) ≈ one_state(1)
    @test apply(zero_state(1),Op("identity transform",identity)) ≈ zero_state(1)
    @test_throws ArgumentError Op("X",1;if1="X")
    @test !isdefined(BlueTangle,:measure_ZNE)
end

@testset "Callback sampling and controlled routing" begin
    calls = Ref(0)
    callback = Op("alternating preparation", state -> begin
        calls[] += 1
        isodd(calls[]) ? apply(state, Op("X",1)) : state
    end)
    operations = [Op("I",1), callback]
    result = measure(operations,10;backend=:density_matrix)
    @test calls[] == 10
    @test result.expect == [0.0]
    op = Op("X",1;control=4)
    routed = route([op], Layout(ones(Int,1,4)))
    @test sa.sparse(routed) ≈ op.expand(4)
    @test all(g -> g isa Op && (g.q == 1 ? abs(g.control-g.qubit)==1 : abs(g.qubit-g.target_qubit)==1), routed)
end

@testset "Measurements inside conditional branches" begin
    conditional = ifOp("MZ", 1, [Op("MX", 2), Op("MZ", 1)], Op("I", 2))
    for outer_bit in (0, 1), seed in 1:4
        vector = outer_bit == 0 ? ComplexF64[1, 0, 0, 0] : ComplexF64[0, 0, 1, 0]
        rho = vector * vector'
        for initial in (vector, sa.sparse(vector), rho, sa.sparse(rho),
                        to_MPS(vector, N_MPS(2); cutoff=0.0))
            Random.seed!(seed)
            actual, bits = apply(initial, conditional; track_measurements=true)
            state = actual isa it.MPS ? to_state(actual) : actual
            actual_rho = state isa AbstractVector ? state * state' : state
            if outer_bit == 0
                # The outer result precedes both measurements in the chosen branch.
                @test bits in ([0, 0, 0], [0, 1, 0])
                if length(bits) == 3
                    x_eigenstate = ComplexF64[1, (-1)^bits[2]] / sqrt(2)
                    expected = kron(ComplexF64[1, 0], x_eigenstate)
                    @test actual_rho ≈ expected * expected' atol=1e-12
                end
            else
                # No measurements from the unselected branch enter the record.
                @test bits == [1]
                @test actual_rho ≈ rho atol=1e-12
            end
        end
    end
    # The direct low-level interface retains its original scalar outcome.
    _, bit = conditional.born_apply(ComplexF64[1, 0, 0, 0])
    @test bit === 0
end

@testset "Reset of an entangled qubit" begin
    vector = ComplexF64[sqrt(0.3), 0, 0, exp(0.4im)*sqrt(0.7)]
    rho = vector * vector'
    identity = Matrix{ComplexF64}(I, 2, 2)
    for qubit in 1:2, name in ("RES", "RESET")
        # Reset has Kraus operators |0><0| and |0><1| on the target qubit.
        local_kraus = (ComplexF64[1 0; 0 0], ComplexF64[0 1; 0 0])
        kraus = [qubit == 1 ? kron(K, identity) : kron(identity, K) for K in local_kraus]
        expected_rho = sum(K * rho * K' for K in kraus)
        for initial in (rho, sa.sparse(rho))
            actual, bits = apply(initial, Op(name, qubit); track_measurements=true)
            @test actual ≈ expected_rho atol=1e-12
            @test isempty(bits)
        end
        for seed in 1:4
            Random.seed!(seed)
            outcome = rand() < 0.3 ? 0 : 1
            expected = zeros(ComplexF64, 4)
            expected[outcome == 0 ? 1 : (qubit == 1 ? 2 : 3)] = 1
            for initial in (vector, sa.sparse(vector), to_MPS(vector, N_MPS(2); cutoff=0.0))
                Random.seed!(seed)
                actual, bits = apply(initial, Op(name, qubit); track_measurements=true)
                state = actual isa it.MPS ? to_state(actual) : actual
                # The target is zero; its partner retains the sampled correlation.
                @test state * state' ≈ expected * expected' atol=1e-12
                @test isempty(bits)
            end
        end
    end
end

@testset "Conditional MPS truncation settings" begin
    branch = [Op("H", 1), Op("CX", 1, 2)]
    conditional = ifOp("MZ", 3, branch, Op("I", 3))
    expected = apply(zero_state(N_MPS(3)), branch; maxdim=1, cutoff=0.0)
    for track in (false, true)
        result = apply(zero_state(N_MPS(3)), conditional;
                       maxdim=1, cutoff=0.0, track_measurements=track)
        actual = track ? first(result) : result
        @test maxlinkdim(actual) == 1
        @test fidelity(to_state(actual), to_state(expected)) ≈ 1
        track && @test last(result) == [0]
    end
end

@testset "Tuple shorthand operations" begin
    bell_ops = [Op("H", 1), Op("CX", 1, 2)]
    bell_tup = [("H", 1), ("CX", 1, 2)]
    bell_mix = [("H", 1), Op("CX", 1, 2)] # Vector{Any}
    bell_ref = apply(bell_ops)

    @test get_N(bell_tup) == 2
    @test get_N(bell_mix) == 2
    @test get_N(Tuple[]) == 0
    @test get_stats(bell_tup) == get_stats(bell_ops)
    @test get_stats(bell_mix) == get_stats(bell_ops)
    @test length(get_layers(bell_tup)) == length(get_layers(bell_ops))
    @test isempty(get_optimized_layers([("X", 1), ("X", 1)]))
    @test to_qasm(bell_tup) == to_qasm(bell_ops)
    @test to_qasm(bell_mix) == to_qasm(bell_ops)
    @test hilbert(bell_tup) ≈ hilbert(bell_ops)
    @test hilbert(bell_mix) ≈ hilbert(bell_ops)

    layout = Layout(ones(Int, 1, 3))
    @test apply(zero_state(3), route([("CX", 1, 3)], layout)) ≈ apply(zero_state(3), Op("CX", 1, 3))

    # Block OpF, single-tuple apply, and * operator across backends
    block = Op("bell_block", bell_tup)
    @test block isa OpF
    @test apply(zero_state(2), block) ≈ bell_ref
    @test apply(("X", 2)) ≈ apply(Op("X", 2))
    @test apply(zero_state(2), ("X", 2)) ≈ apply(zero_state(2), Op("X", 2))
    @test (("X", 2) * zero_state(2; backend=:density)) ≈ (Op("X", 2) * zero_state(2; backend=:density))
    @test to_state(apply(zero_state(2; backend=:camps), ("X", 2))) ≈ apply(Op("X", 2))

    # Mid-circuit measurements & reset via tuple shorthand and mixed lists
    mid_tup = [("X", 1), ("MZ", 1), ("RESET", 1), ("MZ", 1)]
    cond_mix = [("X", 1), ifOp("MZ", 1, "I", "X")]
    for backend in (:statevector, :density_matrix, :tensor, :camps)
        st = apply(bell_tup; backend=backend)
        st_m = apply(bell_mix; backend=backend)
        @test expect(st, "Z") ≈ expect(bell_ref, "Z") atol=1e-12
        @test expect(st_m, "Z") ≈ expect(bell_ref, "Z") atol=1e-12

        m = measure(bell_tup, 60; backend=backend)
        @test all(b -> b in (0, 3), m.bitstr)
        @test isempty(m.mid_measurements)
        @test all(b -> b in (0, 3), measure(bell_mix, 40; backend=backend).bitstr)
        @test all(b -> b in (0, 3), sample(bell_tup, 25; backend=backend))
        @test all(b -> b in (0, 3), sample(bell_mix, 25; backend=backend))

        m_mid = measure(mid_tup, 3; backend=backend)
        @test m_mid.mid_measurements == fill([1, 0], 3)
        @test m_mid.bitstr == [0]

        m_cond = measure(cond_mix, 3; backend=backend)
        @test m_cond.mid_measurements == fill([1], 3)
        @test m_cond.bitstr == [0]

        @test measure(Tuple[], 4; N=2, backend=backend).bitstr == [0]
        @test sample(Tuple[], 4; N=2, backend=backend) == fill(0, 4)
    end

    # CAMPS helper with RESET tuple (returns OpQC, not Op)
    @test to_state(camps([("X", 1), ("RESET", 1)])) ≈ zero_state(1)

    # Shadow with tuple shorthand
    for backend in (:statevector, :tensor, :camps)
        snap = shadow([("I", 1)], 1; backend=backend)
        @test snap isa sa.SparseMatrixCSC
        @test eigvals(Hermitian(Matrix(snap))) ≈ [-1, 2]
    end

    # Text plotting and variational_apply with tuple shorthands
    io = IOBuffer()
    plotqt(bell_tup; io=io)
    @test occursin("H", String(take!(io)))
    plotqt(("X", 1); io=io)
    @test occursin("X", String(take!(io)))
    @test variational_apply([0.4], [("RX", _RX, 1), ("CX", 1, 2)]) ≈
          apply(zero_state(2), [Op("RX(0.4)", 1), Op("CX", 1, 2)])

    # Untyped numeric statevectors still work; invalid Vector{Any} errors cleanly
    @test get_N(Any[1.0, 0.0]) == 1
    @test measure(Any[1.0, 0.0], 5).bitstr == [0]
    @test sample(Any[1.0, 0.0], 5) == fill(0, 5)
    @test_throws ArgumentError apply(Any[1, 2])
    @test_throws ArgumentError measure(Any["not_an_op"], 5)
    @test_throws MethodError measure(bell_tup)
end
