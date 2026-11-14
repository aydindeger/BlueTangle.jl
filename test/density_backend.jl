@testset "Density-matrix evolution" begin
    plus = ComplexF64[1, 1] / sqrt(2)
    plus_y = ComplexF64[1, im] / sqrt(2)
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    rho_plus = plus * plus'
    mixed = Matrix{ComplexF64}(I, 2, 2) / 2
    zero = ComplexF64[1 0; 0 0]
    for storage in (identity, sa.sparse)
        rho = storage(rho_plus)
        original = copy(rho)
        @test apply(rho, Op("X", 1)) ≈ rho
        @test apply(rho, Op("MZ", 1)) ≈ mixed
        @test apply(rho, Op("MX", 1)) ≈ rho
        @test apply(rho, Op("MY", 1)) ≈ mixed
        @test apply(rho, Op("MR", 1)) ≈ (rho_plus + 2mixed) / 3
        @test apply(rho, ifOp("MZ", 1, "I", "X")) ≈ zero
        @test apply(storage(mixed), ifOp("MX", 1, "I", "Z")) ≈ rho_plus
        @test apply(storage(mixed), ifOp("MY", 1, "I", "Z")) ≈ plus_y * plus_y'
        @test apply(rho, ifOp("MR", 1)) ≈ apply(rho, Op("MR", 1))
        @test apply(rho, Op("RESET", 1)) ≈ zero
        @test born_measure_Z(1, rho, 1) ≈ mixed
        @test rho == original

        # A selective measurement samples one branch and normalizes by its trace.
        selected, bits = apply(storage(zero), Op("MZ", 1); track_measurements=true)
        @test selected ≈ zero
        @test bits == [0]
        selected, bits = apply(storage(plus_y * plus_y'), Op("MY", 1); track_measurements=true)
        @test selected ≈ plus_y * plus_y'
        @test bits == [0]
        selected, bits = apply(storage(zero), [Op("X", 1), ifOp("MZ", 1, "I", "X")]; track_measurements=true)
        @test selected ≈ zero
        @test bits == [1]

        @test apply(rho, Op("I", 1); noise=NoiseModel(false, false)) ≈ rho
        @test apply(rho, Op("I", 1); noise=NoiseModel(false, Noise2("depolarizing", 0.2))) ≈ rho
        damping = OpQC("amplitude_damping", 0.3, 1)
        excited = storage(ComplexF64[0 0; 0 1])
        @test apply(excited, damping) ≈ [0.3 0; 0 0.7]
        @test damping.prob(excited) ≈ [0.7, 0.3]
        @test apply(excited, Op("I", 1); noise=NoiseModel(Noise1("amplitude_damping", 0.3), false)) ≈ [0.3 0; 0 0.7]

        rotation = Op("rotation", _RY, 1)
        @test apply(storage(zero), rotation; pars=[pi]) ≈ ComplexF64[0 0; 0 1]
        @test apply(storage(zero), OpF("X", gate.X)) ≈ ComplexF64[0 0; 0 1]
        @test apply(storage(zero), OpF("sequence", [Op("H", 1), Op("MZ", 1)])) ≈ mixed

        rho_bell = storage(bell * bell')
        reset = apply(rho_bell, ifOp("MZ", 1, "I", "X"))
        @test reset ≈ kron(zero, mixed)
        @test tr(reset) ≈ 1
        @test apply(rho_bell, Op("CX", 1, 2); noise=NoiseModel(Noise1("depolarizing", 0.2), false)) ≈ kron(rho_plus, zero)
        @test mag_moments(rho_bell, "Z", 0) ≈ 1
        @test mag_moments(rho_bell, "Z", 2) ≈ 4
        @test mag_moments(rho_bell, "X", 2) ≈ 4
        @test mag_moments(storage(Matrix{ComplexF64}(I, 4, 4) / 4), "Z"; max_moment=2) ≈ [0, 2]
        @test entanglement_entropy(rho_bell)[1] ≈ log(2)
        @test entanglement_entropy(rho_bell; spectrum_bool=false) ≈ log(2)
        @test entanglement_entropy(rho_bell; spectrum_bool=true)[2] ≈ [log(2), log(2)]

        ghz = ComplexF64[1, 0, 0, 0, 0, 0, 0, 1] / sqrt(2)
        @test entanglement_entropy(storage(ghz * ghz'); spectrum_bool=false) ≈ log(2)
    end
    @test entanglement_entropy(bell)[1] ≈ log(2)
    @test entanglement_entropy(bell; spectrum_bool=false) ≈ log(2)

    # An asymmetric two-qubit channel must respect the argument order, including
    # reversed and separated qubits, when calculating its branch probabilities.
    gamma = 0.3
    kraus = [kron([1.0 0; 0 sqrt(1-gamma)], gate.I),
             kron([0 sqrt(gamma); 0 0], gate.I)]
    state = product_state([0, 0, 1])
    channel = OpQC("first-qubit damping", kraus, 3, 1)
    @test channel.prob(state) ≈ [0.7, 0.3]
    @test channel.prob(state * state') ≈ [0.7, 0.3]
    @test apply(state * state', channel) ≈ 0.7state*state' + 0.3zero_state(3)*zero_state(3)'
    identity_channel = OpQC("three qubits", [Matrix{Float64}(I, 8, 8)], 1)
    @test apply(Matrix(state * state'), identity_channel) ≈ state * state'
end

@testset "Completed operator constructors" begin
    identity_noise = custom_noise(1, "identity", [Matrix{Float64}(I, 2, 2)])
    @test identity_noise.q == 1
    @test identity_noise.apply(ComplexF64[1 0; 0 0], 1) ≈ [1 0; 0 0]
    @test_throws ArgumentError custom_noise(2, "wrong size", [Matrix{Float64}(I, 2, 2)])
    @test to_state(apply(zero_state(N_MPS(1)), OpF("X", gate.X))) ≈ one_state(1)
    for op in (Op("RY", _RY, 1), Op("U3", _U3, 1))
        pars = op.name == "RY" ? [0.3] : [0.3, -0.2, 0.4]
        @test adjoint(op).expand(1, pars...) ≈ op.expand(1, pars...)'
    end
    @test_throws ArgumentError Op("CX", 1, 2; control=1)
    @test_throws ArgumentError Op("CX", 1, 2; control=2)

    sites = N_MPS(3)
    state = normalize(ComplexF64[1, 2im, -2, 3, im, 1-im, 0.2, -1])
    mps = to_MPS(state, sites)
    for op in (Op("RY(0.3)", 1; control=3),
               Op("custom", _RX(0.2), 3; control=1),
               Op("CX", 1, 2; control=3),
               Op("SWAP", 2, 3; control=1))
        @test to_state(apply(mps, op)) ≈ apply(state, op)
    end
    op = Op("parameterized", _RY, 1; control=3)
    @test to_state(apply(mps, op; pars=[0.3])) ≈ apply(state, op; pars=[0.3])

    # MR must undo exactly the basis rotation sampled for that measurement.
    # Applying the same randomly chosen measurement twice is idempotent.
    for seed in 1:6
        Random.seed!(seed)
        once = apply(mps, Op("MR", 1))
        Random.seed!(seed)
        twice = apply(once, Op("MR", 1))
        @test fidelity(once, twice) ≈ 1
    end
end

@testset "MPS channels without statevector expansion" begin
    sites = N_MPS(40)
    psi = one_state(sites)
    damped = apply(psi, OpQC("amplitude_damping", 1.0, 20))
    @test expect(damped, "Z")[20] ≈ 1
    @test expect(damped, "Z")[19] ≈ -1
    damped_pair = apply(psi, OpQC("amplitude_damping", 1.0, 40, 1))
    @test expect(damped_pair, "Z")[[1, 40]] ≈ [1, 1]
    three_qubit = OpQC("flip first", [kron(gate.X, gate.I, gate.I)], 4)
    flipped = apply(psi, three_qubit)
    @test expect(flipped, "Z")[[4, 5, 6]] ≈ [1, -1, -1]
    small_sites = N_MPS(2)
    @test BlueTangle.calculate_high_moments(plus_state(small_sites), small_sites, "Z", 2) ≈ [0, 2]
end

@testset "Nonlocal controlled custom gates" begin
    sites = N_MPS(4)
    state = normalize(ComplexF64[complex(i, (-1)^i) for i in 1:16])
    psi = to_MPS(state, sites)
    for op in (Op("SWAP", 1, 4; control=2),
               Op("custom", _FSIM(0.2, 0.3), 4, 1; control=3),
               Op("parameterized", _FSIM, 1, 4; control=3))
        pars = op.type == "f" ? (; pars=[0.2, 0.3]) : (;)
        expected = apply(state, op; pars...)
        @test to_state(apply(psi, op; pars...)) ≈ expected
        @test norm(expected) ≈ 1
    end
end

@testset "Data persistence" begin
    operations = [Op("H", 1)]
    options = Options(circuit_name="saved operations")
    measurement = measure(plus_state(1))
    mktempdir() do directory
        path = joinpath(directory, "result.ser")
        @test BlueTangle.save_data(operations, measurement, ["note"]; options=options, path=path) == path
        loaded = BlueTangle.load_data(path)
        @test loaded.other == ["note"]
        @test loaded.measurement.sample ≈ measurement.sample
        @test apply(zero_state(1), loaded.ops) ≈ plus_state(1)
        @test loaded.options.circuit_name == options.circuit_name
        @test BlueTangle.load_data(path[1:end-4]).other == ["note"]
    end
end

@testset "Channel names preserve supplied Kraus operators" begin
    excited = one_state(1)
    for name in ("reset", "RES")
        channel = custom_noise(1, name, [Matrix{Float64}(I, 2, 2)])
        @test_throws ArgumentError to_qasm([Op(channel, 1)])
        for initial in (excited, to_rho(excited), to_MPS(excited, N_MPS(1)))
            result, bits = apply(initial, Op(channel, 1); track_measurements=true)
            @test isempty(bits)
            @test initial isa it.MPS ? to_state(result) ≈ excited : result ≈ initial
        end
    end

    # Built-in reset still prepares |0>, without recording a measurement bit.
    for name in ("RESET", "RES")
        @test occursin("reset q[0];", to_qasm([Op(name, 1)]))
        for initial in (excited, to_rho(excited), to_MPS(excited, N_MPS(1)))
            result, bits = apply(initial, [Op(name, 1), Op("MZ", 1)]; track_measurements=true)
            @test bits == [0]
            @test initial isa AbstractMatrix ? result ≈ to_rho(zero_state(1)) :
                  (initial isa it.MPS ? to_state(result) : result) ≈ zero_state(1)
        end
    end

    # Reset on an entangled state retains the untouched qubit's mixed ensemble.
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    expected = kron(ComplexF64[1 0; 0 0], Matrix{ComplexF64}(I, 2, 2) / 2)
    for initial in (bell, to_rho(bell), to_MPS(bell, N_MPS(2)))
        result, bits = apply(initial, Op("RESET", 1); track_measurements=true)
        @test isempty(bits)
        if initial isa AbstractMatrix
            @test result ≈ expected
        else
            vector = initial isa it.MPS ? to_state(result) : result
            @test norm(vector) ≈ 1
            @test max(abs2(vector[1]), abs2(vector[2])) ≈ 1
        end
    end
end
