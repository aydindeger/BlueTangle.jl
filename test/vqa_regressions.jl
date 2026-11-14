@testset "Training edge cases and gradient semantics" begin
    given32 = Float32[0.31]
    rotation = QuantumOps[Op("RY", _RY, 1)]
    energy(p) = real(expect(variational_apply(p, rotation), gate.Z))
    @test only(BlueTangle.ForwardDiff.gradient(energy, given32)) ≈ -sin(only(given32)) atol=2e-6
    result32 = train(energy, given32; optimizer=:descent, iterations=2, history=false)
    @test result32.loss ≈ cos(only(result32.parameters)) atol=2e-6
    @test given32 == Float32[0.31]

    # A fixed circuit needs one objective evaluation, regardless of optimizer.
    for optimizer in (:adam, :descent, :radam, :momentum, :nesterov, :lbfgs,
                      :prima, :cobyla, :uobyqa, :newuoa, :bobyqa, :lincoa)
        calls = Ref(0)
        fixed_loss(p) = begin
            calls[] += 1
            real(expect(variational_apply(p, [Op("H", 1)]), gate.Z))
        end
        result = train(fixed_loss, Float64[]; optimizer=optimizer, iterations=2, history=false)
        @test result.loss ≈ 0.0 atol=1e-12
        @test isempty(result.parameters)
        @test calls[] == 1
    end

    # Each optimization starts with fresh optimizer state and leaves its inputs
    # unchanged, even for algorithms with momentum or running moments.
    for optimizer in (:adam, :descent, :radam, :momentum, :nesterov)
        initial_parameters = [0.31]
        initial_state = ComplexF64[1, 0]
        objective(p) = real(expect(variational_apply(p, rotation; init=initial_state), gate.Z))
        first_result = train(objective, initial_parameters; optimizer=optimizer, iterations=3)
        second_result = train(objective, initial_parameters; optimizer=optimizer, iterations=3)
        @test first_result == second_result
        @test last(first_result.history) ≈ cos(only(first_result.parameters))
        @test initial_parameters == [0.31]
        @test initial_state == ComplexF64[1, 0]
    end

    # Differentiating a sampled Kraus branch omits the derivative of its Born
    # probability. The exact damping-channel expectation is gamma+(1-gamma)cosθ.
    theta = [1.1]
    gamma = 0.4
    noise = NoiseModel(Noise1("amplitude_damping", gamma), false)
    rho0 = ComplexF64[1 0; 0 0]
    density_loss(p) = expect(variational_apply(p, rotation; init=rho0, noise=noise), gate.Z)
    @test density_loss(theta) ≈ gamma + (1-gamma)*cos(only(theta))
    @test only(BlueTangle.ForwardDiff.gradient(density_loss, theta)) ≈ -(1-gamma)*sin(only(theta))
    @test rho0 == ComplexF64[1 0; 0 0]

    density_result = train(density_loss, theta; optimizer=:descent, iterations=2, history=false)
    @test density_result.loss ≈ gamma + (1-gamma)*cos(only(density_result.parameters))

    for initial_state in (zero_state(1), zero_state(N_MPS(1)))
        trajectory_loss(p) = expect(variational_apply(p, rotation; init=initial_state, noise=noise), Op("Z", 1))
        @test_throws ArgumentError train(trajectory_loss, theta; optimizer=:descent, iterations=1)
    end

    for stochastic_op in (Op("MZ", 1), OpQC("amplitude_damping", gamma, 1),
                          ifOp("MZ", 1, Op("I", 1), Op("X", 1)))
        ops = QuantumOps[rotation[1], stochastic_op]
        sampled_loss(p) = real(expect(variational_apply(p, ops), gate.Z))
        @test_throws ArgumentError train(sampled_loss, theta; optimizer=:descent, iterations=1)
    end
end

@testset "Variational conditional noise" begin
    ops = QuantumOps[ifOp("MZ", 1, Op("I", 2), Op("X", 2))]
    noise = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    initial_state = product_state([1, 0])
    # The selected X gate is followed by damping, as in ordinary application.
    @test variational_apply(Float64[], ops; init=initial_state, noise=noise) ≈ apply(initial_state, ops; noise=noise)
    @test initial_state == product_state([1, 0])
end
