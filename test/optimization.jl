@testset "Twirling preserves operated qubits" begin
    function unitary_of(ops, N)
        U = Matrix{ComplexF64}(I, 2^N, 2^N)
        for op in ops
            U = op.expand(N) * U
        end
        return U
    end

    Random.seed!(1307)
    for name in ("CX", "CNOT", "ECR"), (q, t) in ((1, 3), (3, 1))
        op = Op(name, q, t)
        expected = Matrix(op.expand(3))
        for _ in 1:32
            twirled = apply_twirl(op)
            result = unitary_of(twirled, 3)
            phase = tr(expected' * result) / 8
            @test abs(phase) ≈ 1
            @test result ≈ phase * expected
            @test all(g -> g === op || g.qubit in (q, t), twirled)
        end
    end

    controlled = Op("CX", 1, 2; control=3)
    @test only(apply_twirl(controlled)) === controlled
    @test !BlueTangle._is_cnot_like(controlled)
    custom = Op("CX", Matrix{ComplexF64}(I, 4, 4), 1, 2)
    @test only(apply_twirl(custom)) === custom
    @test !BlueTangle._is_cnot_like(custom)

    # The common sampling path must also leave Toffoli gates intact.
    operations = [Op("H", 1), Op("X", 2), controlled]
    expected = apply(zero_state(3), operations)
    result = measure(operations, 1; options=Options(twirl=true, density_matrix=true),
                     backend=:density_matrix)
    @test result.density_matrix ≈ expected * expected'
end

@testset "Conservative ideal gate cancellation" begin
    for name in ("X", "Y", "Z", "H")
        @test isempty(get_optimized_layers([Op(name, 1), Op(name, 1)]))
    end
    for name in ("S", "SX", "XSQRT")
        @test isempty(get_optimized_layers([Op(name, 1) for _ in 1:4]))
    end
    for name in ("CX", "CNOT", "CY", "CZ", "ECR", "SWAP", "FSWAP")
        @test isempty(get_optimized_layers([Op(name, 1, 2), Op(name, 1, 2)]))
    end

    conditional = ifOp("MZ", 1)
    function_op = OpF("function", identity)
    @test only(only(get_optimized_layers([conditional]))) === conditional
    @test only(only(get_optimized_layers([function_op]))) === function_op
    @test length(get_optimized_layers([Op("I", 1), function_op])) == 2

    # Labels are user-controlled; a label alone cannot establish an identity.
    custom = [Op("X", gate.X, 1), Op("X", gate.H, 1)]
    optimized = get_optimized_layers(custom)
    @test length(optimized) == 2
    @test apply(vcat(optimized...), zero_state(1)) ≈ minus_state(1)
    @test length(get_optimized_layers([Op("X", gate.H, 1), Op("X", gate.X, 1)])) == 2
    @test length(get_optimized_layers([Op("CX", gate.CX, 1, 2), Op("CX", gate.SWAP, 1, 2)])) == 2

    parameterized = Op("X", θ -> gates("RX($θ)"), 1)
    @test length(get_optimized_layers([parameterized, parameterized])) == 2
    @test length(get_optimized_layers([Op("X", 1; noisy=false), Op("X", 1; noisy=true)])) == 2
    @test isempty(get_optimized_layers([Op("X", 1; noisy=false), Op("X", 1; noisy=false)]))
end
