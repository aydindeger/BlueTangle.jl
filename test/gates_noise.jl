@testset "Gate construction and Pauli noise" begin
    theta = 0.37
    for (name, pauli) in (("RX", gate.X), ("RY", gate.Y), ("RZ", gate.Z))
        matrix, arguments, values = gates("$name($theta)", true)
        @test matrix ≈ exp(-im * theta / 2 * pauli)
        @test arguments == "($theta)"
        @test values == [theta]
    end
    @test gates("RX(pi/2)") ≈ exp(-im * pi / 4 * gate.X)
    @test gates("h") == gate.H
    @test gates("FSIM(0.0,0.0)") == Matrix{Float64}(I, 4, 4)
    @test gates("U3(0.0,0.0,0.0)", true)[2:3] == ("(0.0,0.0,0.0)", [0.0, 0.0, 0.0])

    state = normalize(ComplexF64[1, 2im])
    rho = state * state'
    for (model, pauli) in (("phase_flip", gate.Z), ("bit_flip", gate.X), ("bit_phase_flip", gate.Y))
        for p in (0.0, 0.23, 1.0)
            kraus = noise_model1(model, p)
            @test is_valid_quantum_channel(kraus)
            @test sum(k * rho * k' for k in kraus) ≈ (1 - p) * rho + p * pauli * rho * pauli'
        end
    end
    @test noise_model1("PHASE_FLIP", 1.0)[2] == gate.Z
    @test noise_model2("phase_flip", 1.0)[end] == kron(gate.Z, gate.Z)
    @test apply(zero_state(1), OpQC("phase_flip", 1.0, 1)) ≈ zero_state(1)
end

@testset "Custom channel validation and Kraus ownership" begin
    identity2 = Matrix{ComplexF64}(I, 2, 2)
    for kraus in (Matrix{Float64}[], [zeros(0, 0)], [ones(2, 3)],
                  [ones(2, 2), ones(4, 4)], [fill(NaN, 2, 2)],
                  [fill(Inf, 2, 2)], [0.5identity2])
        @test !is_valid_quantum_channel(kraus)
        @test_throws ArgumentError QuantumChannel("invalid", kraus, 0.0)
        @test_throws ArgumentError OpQC("invalid", kraus, 1)
        @test_throws ArgumentError custom_noise(1, "invalid", kraus)
    end
    @test !is_valid_quantum_channel([1, 2])
    @test !is_valid_quantum_channel([fill("invalid", 2, 2)])
    @test is_valid_quantum_channel([identity2])
    @test is_valid_quantum_channel([sa.sparse(identity2)])
    @test is_valid_quantum_channel([sqrt(0.7)*gate.I, sqrt(0.3)*gate.Z])
    for q in (-1, 0, 3)
        @test_throws ArgumentError QuantumChannel(q, "bit_flip", 0.1)
    end
    for dimension in (1, 3, 16)
        kraus = [Matrix{Float64}(I, dimension, dimension)]
        @test_throws ArgumentError QuantumChannel("unsupported size", kraus, 0.0)
        @test_throws ArgumentError OpQC("unsupported size", kraus, 1)
    end
    @test QuantumChannel(2, "bit_flip", 0.1).q == 2
    @test OpQC("three qubits", [Matrix{Float64}(I, 8, 8)], 1).q == 3

    for construct in (ks -> custom_noise(1, "custom", ks), ks -> OpQC("custom", ks, 1))
        input = [zeros(ComplexF64, 2, 2), copy(identity2)]
        channel = construct(input)
        # Neither replacing an input entry nor editing its matrix changes the
        # constructed channel. All backends use the same owned Kraus operators.
        input[1] = Matrix{ComplexF64}(gate.X)
        input[2] .= 0
        @test channel.kraus[1] == zeros(2, 2)
        @test channel.kraus[2] == identity2

        # Replace the stored identity channel with a deterministic bit flip.
        # Its branch probabilities change too, exposing stale closure captures.
        channel.kraus[1] = Matrix{ComplexF64}(gate.X)
        channel.kraus[2] = zeros(ComplexF64, 2, 2)
        initial = zero_state(1)
        rho = to_rho(initial)
        args = channel isa QuantumChannel ? (1,) : ()
        @test channel.prob(initial, args...) ≈ [1, 0]
        @test channel.prob(rho, args...) ≈ [1, 0]
        @test channel.apply(initial, args...) ≈ one_state(1)
        @test channel.apply(rho, args...) ≈ to_rho(one_state(1))
        op = channel isa QuantumChannel ? Op(channel, 1) : channel
        @test apply(initial, op) ≈ one_state(1)
        @test apply(rho, op) ≈ to_rho(one_state(1))
        @test to_state(apply(zero_state(N_MPS(1)), op)) ≈ one_state(1)
        if channel isa QuantumChannel
            noise = NoiseModel(channel, false)
            @test apply(initial, Op("I", 1); noise) ≈ one_state(1)
            @test apply(rho, Op("I", 1); noise) ≈ to_rho(one_state(1))
            @test to_state(apply(zero_state(N_MPS(1)), Op("I", 1); noise)) ≈ one_state(1)
        end
    end
end
