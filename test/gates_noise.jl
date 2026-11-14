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
