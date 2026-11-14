@testset "Numerical utilities" begin
    x = 1e9 .+ collect(0.0:4.0)
    a, b, se_a, se_b = linear_fit(x, 3 .+ 2 .* x)
    @test a ≈ 3
    @test b ≈ 2
    @test se_a == se_b == 0
    @test isnan(linear_fit([0.0, 1.0], [2.0, 4.0])[3])
    @test_throws ArgumentError linear_fit([1.0, 1.0], [0.0, 1.0])
    x = 1e3 .+ collect(0.0:5.0)
    coefficients = quadratic_fit(x, 2x.^2 + 3x .+ 4)
    @test size(coefficients) == (3, 1)
    @test coefficients ≈ reshape([2, 3, 4], 3, 1) atol=1e-3
    @test same_up_to_global_phase(ComplexF64[2, 0], ComplexF64[2im, 0])[1]
    @test !same_up_to_global_phase(ComplexF64[2, 0], ComplexF64[1im, 0])[1]

    rng = MersenneTwister(19)
    for N in 1:3
        matrix = randn(rng, ComplexF64, 2^N, 2^N)
        coefficients, nonzero, names = pauli_decomposition(matrix)
        @test pauli_reconstruction(coefficients, pauli_decomposition_tensor(N)) ≈ matrix
        @test pauli_reconstruction(nonzero, names) ≈ matrix
    end
    @test pauli_reconstruction(zeros(4), pauli_decomposition_tensor(1)) == zeros(2, 2)
    for _ in 1:10
        A, B = randn(rng, ComplexF64, 2, 2), randn(rng, ComplexF64, 2, 2)
        C = kron(A, B)
        a, b = kronecker_decomposition(C)
        @test kron(a, b) ≈ C
        U = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        alpha, beta, gamma, delta = zyz_decomposition(U)
        @test exp(im * alpha) * _RZ(beta) * _RY(gamma) * _RZ(delta) ≈ U
    end
end

@testset "Relabeling preserves operation data" begin
    custom = Op("custom", _RY(0.3), 1; noisy=false)
    relabeled, mapping = relabel_swap([Op("SWAP", 1, 2), custom])
    @test mapping == [2, 1]
    @test relabeled[1].qubit == 2
    @test relabeled[1].mat == custom.mat
    @test !relabeled[1].noisy
    variable = Op("rotation", _RY, 1)
    relabeled, _ = relabel_swap([Op("SWAP", 1, 2), variable])
    @test relabeled[1].mat === variable.mat
    @test relabeled[1].expand(2, 0.2) ≈ Op("RY(0.2)", 2).expand(2)

    channel = OpQC("custom damping", [[1.0 0; 0 sqrt(0.7)], [0 sqrt(0.3); 0 0]], 1)
    conditional = ifOp("MZ", 1, Op("I", 1), Op("X", 2))
    relabeled, _ = relabel_swap(QuantumOps[Op("SWAP", 1, 2), channel, conditional])
    @test relabeled[1] isa OpQC
    @test relabeled[1].qubit == 2
    @test relabeled[1].kraus == channel.kraus
    @test relabeled[2] isa ifOp
    @test relabeled[2].qubit == 2
    @test relabeled[2].if01[2][1].qubit == 1
    opaque = OpF("opaque", identity)
    @test_throws ArgumentError relabel_swap(QuantumOps[Op("SWAP", 1, 2), opaque])
    @test_throws ArgumentError relabel_swap([Op("SWAP", 1, 2; control=3)])
    fake_swap = Op("SWAP", kron(gate.X, gate.I), 1, 2)
    @test length(first(relabel_swap([fake_swap]))) == 1

    # Moving the absorbed SWAP to the end gives the same physical evolution.
    state = normalize(ComplexF64[1, 2im, -3, 1-im])
    original = QuantumOps[Op("SWAP", 1, 2), custom, channel]
    relabeled, _ = relabel_swap(original)
    @test apply(apply(state * state', relabeled), Op("SWAP", 1, 2)) ≈ apply(state * state', original)
end

@testset "Joint two-qubit Born measurement" begin
    for bits in ([0, 0], [1, 0], [0, 1], [1, 1])
        state = product_state(bits)
        result, outcome = born_measure_Z(2, state, 1, 2)
        @test result ≈ state
        @test outcome == 1 + bits[1] + 2bits[2]
    end
end
