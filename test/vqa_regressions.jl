@testset "VQA edge cases and gradient semantics" begin
    # A shift rule for an expectation does not differentiate an arbitrary
    # nonlinear function of that expectation: cos(theta)^2 is a counterexample.
    @test_throws ArgumentError begin
        opt = AnsatzOptions(N=1, ops=["RY"], loss=s -> expect(s, gate.Z)^2,
                            pars_initial=[0.31], paramshift=true,
                            model="descent", number_of_iterations=1)
        VQA(opt)
    end

    given32 = Float32[0.31]
    shift_opt = AnsatzOptions(N=1, ops=["RY"], loss=gate.Z,
                              pars_initial=given32, paramshift=true,
                              model="descent", number_of_iterations=2, history=false)
    _, shift_gradient = BlueTangle.loss_and_grad_paramshift(given32, shift_opt)
    @test only(shift_gradient) ≈ -sin(only(given32)) atol=2e-6
    energy32, parameters32 = VQA(shift_opt)
    @test energy32 ≈ cos(only(parameters32)) atol=2e-6
    @test given32 == Float32[0.31]

    # A fixed circuit needs one objective evaluation, regardless of optimizer.
    for model in ("adam", "descent", "radam", "momentum", "nesterov", "lbfgs",
                  "prima", "cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa")
        opt = AnsatzOptions(N=1, ops=["H"], loss=gate.Z, model=model,
                            number_of_iterations=2, history=false)
        energy, parameters = VQA(opt)
        @test energy ≈ 0.0 atol=1e-12
        @test isempty(parameters)
    end

    # Each optimization starts with fresh optimizer state and leaves its inputs
    # unchanged, even for algorithms with momentum or running moments.
    for model in ("adam", "descent", "radam", "momentum", "nesterov")
        initial_parameters = [0.31]
        initial_state = ComplexF64[1, 0]
        opt = AnsatzOptions(N=1, ops=["RY"], loss=gate.Z, model=model,
                            init=initial_state, pars_initial=initial_parameters,
                            number_of_iterations=3, history=true)
        first_result = VQA(opt)
        second_result = VQA(opt)
        @test first_result == second_result
        energies, parameters, parameter_history = first_result
        @test last(energies) ≈ cos(only(parameters))
        @test last(parameter_history) == parameters
        @test initial_parameters == opt.pars_initial == [0.31]
        @test initial_state == opt.init == ComplexF64[1, 0]
    end

    # Differentiating a sampled Kraus branch omits the derivative of its Born
    # probability. The exact damping-channel expectation is gamma+(1-gamma)cosθ.
    theta = [1.1]
    gamma = 0.4
    noise = NoiseModel(Noise1("amplitude_damping", gamma), false)
    rotation = QuantumOps[Op("RY", _RY, 1)]
    rho0 = ComplexF64[1 0; 0 0]
    density_loss(p) = expect(variational_apply(p, rotation; init=rho0, noise=noise), gate.Z)
    @test density_loss(theta) ≈ gamma + (1-gamma)*cos(only(theta))
    @test only(BlueTangle.ForwardDiff.gradient(density_loss, theta)) ≈ -(1-gamma)*sin(only(theta))
    @test rho0 == ComplexF64[1 0; 0 0]

    density_opt = AnsatzOptions(N=1, ops=rotation, loss=gate.Z, init=rho0,
                                noise=noise, pars_initial=theta, model="descent",
                                number_of_iterations=2, history=false)
    density_energy, density_parameters = VQA(density_opt)
    @test density_energy ≈ gamma + (1-gamma)*cos(only(density_parameters))

    for initial_state in (zero_state(1), zero_state(N_MPS(1)))
        trajectory_opt = AnsatzOptions(N=1, ops=rotation, loss=s -> expect(s, Op("Z", 1)),
                                       init=initial_state, noise=noise, pars_initial=theta,
                                       model="descent", number_of_iterations=1, history=false)
        @test_throws ArgumentError VQA(trajectory_opt)
    end

    for stochastic_op in (Op("MZ", 1), OpQC("amplitude_damping", gamma, 1),
                          ifOp("MZ", 1, Op("I", 1), Op("X", 1)))
        opt = AnsatzOptions(N=1, ops=QuantumOps[rotation[1], stochastic_op],
                            loss=gate.Z, pars_initial=theta, model="descent",
                            number_of_iterations=1, history=false)
        @test_throws ArgumentError VQA(opt)
    end
end

@testset "Variational conditional noise" begin
    ops = QuantumOps[ifOp("MZ", 1, Op("I", 2), Op("X", 2))]
    noise = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    initial_state = product_state([1, 0])
    opt = AnsatzOptions(N=2, ops=ops, loss=hilbert(2, gate.Z, 2),
                        init=initial_state, noise=noise)
    # The selected X gate is followed by damping, as in ordinary application.
    @test variational_apply(Float64[], opt) ≈ apply(ops, initial_state; noise=noise)
    @test initial_state == product_state([1, 0])
end
