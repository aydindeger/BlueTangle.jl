@testset "Plain training API" begin
    # A nonlinear quantum-regression loss has a different derivative from the
    # expectation itself. Check one update against its analytic derivative.
    operations = QuantumOps[Op("RY", _RY, 1), Op("RZ", _RZ, 1)]
    initial = [0.31, -0.4]
    target = 0.2
    objective(p) = (expect(variational_apply(p, operations), gate.Z)-target)^2 + 0.1p[2]^4
    analytic_gradient(p) = [-2(cos(p[1])-target)*sin(p[1]), 0.4p[2]^3]
    @test BlueTangle.ForwardDiff.gradient(objective, initial) ≈ analytic_gradient(initial)
    result = train(objective, initial; optimizer=:descent, iterations=1, learning_rate=0.05)
    @test keys(result) == (:loss, :parameters, :history)
    @test result.parameters ≈ initial - 0.05analytic_gradient(initial)
    @test result.loss ≈ (cos(result.parameters[1])-target)^2 + 0.1result.parameters[2]^4
    @test result.history == [result.loss]
    @test initial == [0.31, -0.4]
    @test result.parameters !== initial
    analytic_result = train(objective, initial; optimizer=:descent, iterations=1,
                            learning_rate=0.05, gradient=analytic_gradient)
    @test analytic_result.parameters ≈ result.parameters

    callback_state = ComplexF64[1, 0]
    callback_ops = QuantumOps[Op("RY", _RY, 1), Op("Hadamard callback", s -> gate.H*s)]
    callback_loss(p) = expect(variational_apply(p, callback_ops; init=callback_state), gate.Z)
    @test only(BlueTangle.ForwardDiff.gradient(callback_loss, [0.31])) ≈ cos(0.31)
    callback_result = train(callback_loss, [0.31]; optimizer=:descent, iterations=1)
    @test only(callback_result.parameters) ≈ 0.31 - 0.01cos(0.31)
    @test callback_state == ComplexF64[1, 0]

    # Every optimizer exposes the same fields, including when history is off.
    # PRIMA records trial evaluations; its last trial need not be its optimum.
    quadratic(p) = (p[1]-0.2)^2 + (p[2]+0.7)^2
    for optimizer in (:adam, :descent, :gradient, :radam, :momentum, :nesterov,
                      :lbfgs, :prima, :cobyla, :uobyqa, :newuoa, :bobyqa, :lincoa)
        parameters = [0.8, -0.1]
        fitted = train(quadratic, parameters; optimizer=optimizer, iterations=12)
        @test keys(fitted) == (:loss, :parameters, :history)
        @test fitted.loss ≈ quadratic(fitted.parameters)
        @test fitted.loss <= quadratic(parameters) + 1e-10
        @test !isempty(fitted.history)
        quiet = train(quadratic, parameters; optimizer=optimizer, iterations=12, history=false)
        @test keys(quiet) == (:loss, :parameters, :history)
        @test isempty(quiet.history)
        @test quiet.loss ≈ quadratic(quiet.parameters)
        @test quiet.parameters ≈ fitted.parameters
        @test parameters == [0.8, -0.1]
    end

    calls = Ref(0)
    fixed_loss(p) = (calls[] += 1; sum(abs2, p) + 3.0)
    empty_result = train(fixed_loss, Float64[]; optimizer=:prima, iterations=20)
    @test calls[] == 1
    @test empty_result == (loss=3.0, parameters=Float64[], history=[3.0])
    stopped = train(fixed_loss, [0.5]; iterations=0, history=false)
    @test calls[] == 2
    @test stopped == (loss=3.25, parameters=[0.5], history=Float64[])
    stopped_with_history = train(fixed_loss, [0.5]; iterations=0)
    @test calls[] == 3
    @test stopped_with_history.history == [3.25]
end

@testset "Direct variational states and noise" begin
    operations = QuantumOps[Op("RY", _RY, 1)]
    prepared = apply(zero_state(2), [Op("X", 2)])
    parameter = [0.43]
    result = variational_apply(parameter, operations; init=prepared)
    @test result ≈ kron(_RY(parameter[1])*ComplexF64[1, 0], ComplexF64[0, 1])
    @test prepared == product_state([0, 1])
    @test parameter == [0.43]

    damping = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    fixed = [Op("X", 1)]
    @test variational_apply(Float64[], fixed; noise=damping) ≈ zero_state(1)
    @test variational_apply(Float64[], fixed) ≈ one_state(1)

    conditional = QuantumOps[ifOp("MZ", 1, Op("I", 2), Op("X", 2))]
    initial = product_state([1, 0])
    @test variational_apply(Float64[], conditional; init=initial, noise=damping) ≈ initial
    @test variational_apply(Float64[], conditional; init=initial) ≈ product_state([1, 1])
    @test initial == product_state([1, 0])
end

@testset "MPS training differentiation" begin
    sites = N_MPS(2)
    initial_state = zero_state(sites)
    operations = QuantumOps[Op("RY", _RY, 1), Op("RY", _RY, 2), Op("CX", 1, 2)]
    initial_parameters = [1.1, 0.8]
    # Positive cutoff exercises the tensor decomposition path. Compare the
    # finite-difference update to an independent, smaller central difference.
    objective(p) = expect(variational_apply(p, operations; init=initial_state, cutoff=0.5), Op("Z", 1))
    delta = 1e-6
    reference_gradient = [begin
        offset = zeros(2); offset[i] = delta
        (objective(initial_parameters+offset)-objective(initial_parameters-offset))/(2delta)
    end for i in 1:2]
    result = train(objective, initial_parameters; optimizer=:descent, iterations=1,
                   learning_rate=0.01, gradient=:finite_difference)
    @test result.parameters ≈ initial_parameters - 0.01reference_gradient atol=1e-8
    @test result.loss ≈ objective(result.parameters)
    @test to_state(initial_state) ≈ zero_state(2)
    @test initial_parameters == [1.1, 0.8]

    error = try
        train(objective, initial_parameters; optimizer=:descent, iterations=1)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    @test occursin("finite_difference", sprint(showerror, error))
end

@testset "Ansatz layer construction" begin
    trainable = generate_ansatz_circuit(3, 1, ["RY"])
    fixed = generate_ansatz_circuit(3, 1, ["RY(0.3)"])
    @test all(op -> op.mat ≈ _RY(0.3), fixed)
    expected = kron(_RY(0.3)*[1.0, 0.0], _RY(0.3)*[1.0, 0.0], _RY(0.3)*[1.0, 0.0])
    @test variational_apply(Float64[], fixed) ≈ expected
    @test variational_apply(fill(0.3, 3), trainable) ≈ expected

    deep = generate_ansatz_circuit(5, 2, ["CX"]; deep_circuit=true)
    brick = generate_ansatz_circuit(5, 2, ["CX"]; deep_circuit=false)
    @test [(op.qubit, op.target_qubit) for op in deep] ==
          [(1, 2), (2, 3), (3, 4), (4, 5), (1, 2), (2, 3), (3, 4), (4, 5)]
    @test [(op.qubit, op.target_qubit) for op in brick] == [(1, 2), (3, 4), (2, 3), (4, 5)]
    rotations = EfficientSU2(3, 0)
    @test [(op.name, op.qubit) for op in rotations] ==
          [("RY", 1), ("RY", 2), ("RY", 3), ("RZ", 1), ("RZ", 2), ("RZ", 3)]
    @test isempty(generate_ansatz_circuit(3, 0))
    @test variational_apply(Float64[], QuantumOps[]; N=2) ≈ zero_state(2)

    for constructor in (EfficientSU2, generate_ansatz_circuit)
        @test_throws ArgumentError constructor(0)
        @test_throws ArgumentError constructor(-1)
        @test_throws ArgumentError constructor(2, -1)
    end
end

@testset "Materialized parameter matrices reject stochastic AD" begin
    damping = NoiseModel(Noise1("amplitude_damping", 0.4), false)
    # Parameters can enter through matrices built in the loss closure, without
    # being arguments of function-valued operations.
    objective(p) = expect(variational_apply(Float64[],
        QuantumOps[Op("parameterized matrix", _RY(p[1]), 1)]; noise=damping), gate.Z)
    @test_throws ArgumentError train(objective, [1.1]; optimizer=:descent, iterations=1)

    captured_loss(p) = expect(variational_apply(Float64[],
        QuantumOps[Op("captured parameter", () -> _RY(p[1]), 1)]; noise=damping), gate.Z)
    @test_throws ArgumentError train(captured_loss, [1.1]; optimizer=:descent, iterations=1)

    density_loss(p) = expect(variational_apply(Float64[],
        QuantumOps[Op("parameterized matrix", _RY(p[1]), 1)];
        init=ComplexF64[1 0; 0 0], noise=damping), gate.Z)
    @test only(BlueTangle.ForwardDiff.gradient(density_loss, [1.1])) ≈ -0.6sin(1.1)
end
