@testset "Variational setup and gradients" begin
    given = [0.31]
    ops = [Op("RY", _RY, 1)]
    energy(p) = real(expect(variational_apply(p, ops), gate.Z))
    @test variational_apply(given, ops) ≈ _RY(given[1]) * zero_state(1)
    result = train(energy, given; optimizer=:descent, iterations=3, history=false)
    @test result.loss ≈ cos(only(result.parameters))
    @test result.loss < cos(only(given))
    @test given == [0.31]
    @test result.parameters !== given
    @test_throws DimensionMismatch variational_apply([0.1, 0.2], ops)
    @test_throws DimensionMismatch variational_apply(Float64[], ops)
    @test_throws DimensionMismatch variational_apply([0.1], [Op("U3", _U3, 1)])
    @test_throws ArgumentError train(energy, given; optimizer=:typo)

    prepared = apply(zero_state(1), [Op("X",1)])
    @test variational_apply([0.0], ops; init=prepared) ≈ one_state(1)
    swap_ops = generate_ansatz_circuit(2, 1, ["SWAPA"])
    @test swap_ops[1].q == 2
    @test variational_apply([0.3], swap_ops; init=product_state([0, 1])) ≈
          swap_ops[1].expand(2, 0.3) * product_state([0, 1])
    @test_throws ArgumentError hamiltonian(2,[1.0,"X",2.0])
    @test_throws ArgumentError hamiltonian(2,[1.0,"X"],"typo")
    @test_throws ArgumentError hamiltonian((0,2),[1.0,"X"])
    adjoint_ops = [adjoint(Op("U3", _U3, 1))]
    @test variational_apply([0.2, 0.3, 0.4], adjoint_ops) ≈ _U3(0.2,0.3,0.4)' * zero_state(1)
    adjoint_energy(p) = real(expect(variational_apply(p, adjoint_ops), gate.Z))
    @test BlueTangle.ForwardDiff.gradient(adjoint_energy, [0.2, 0.3, 0.4]) ≈ [-sin(0.2),0,0] atol=1e-10

    # Independent finite differences check every built-in gate's gradient,
    # including gates with a different angle scale and controlled rotations.
    for (name, fn) in BlueTangle._parameterized_gates
        count = BlueTangle._find_argument_number(fn)
        p = collect(range(0.21, step=0.17, length=count))
        two_qubit = size(fn(p...),1) == 4
        for controlled in (false,true)
            n = two_qubit ? 3 : 2
            op = two_qubit ? Op(name,fn,1,2;control=controlled ? 3 : -2) :
                             Op(name,fn,1;control=controlled ? 2 : -2)
            init = normalize(ComplexF64[complex(i,(-1)^i) for i in 1:2^n])
            observable = hilbert(n,gate.X,1) + 0.3hilbert(n,gate.Z,n)
            objective(p) = real(expect(variational_apply(p, [op]; init=init), observable))
            gradient = BlueTangle.ForwardDiff.gradient(objective, p)
            finite_difference = [begin
                step = zeros(length(p)); step[i] = 1e-6
                (objective(p+step)-objective(p-step))/2e-6
            end for i in eachindex(p)]
            @test gradient ≈ finite_difference atol=1e-8
        end
    end

    for optimizer in (:adam, :lbfgs)
        result = train(energy, [0.3]; optimizer=optimizer,
                       iterations=optimizer == :lbfgs ? 12 : 3,
                       history=false, gradient=p -> [-sin(only(p))])
        @test result.loss ≈ cos(only(result.parameters))
    end
end
