@testset "Classical shadows and local reconstruction" begin
    operations = [Op("H", 1), Op("S", 1), Op("CX", 1, 2),
                  Op("RY(0.4)", 3), Op("CX", 2, 3)]

    # Identical measurement records must give the partial trace of the full
    # estimator, including for nonadjacent qubits and an unsorted selection.
    Random.seed!(8101)
    full = shadow(operations, 47)
    @test full isa sa.SparseMatrixCSC
    @test size(full) == (8, 8)
    for qubits in ([1, 3], [3, 1], [2], 2, Int[])
        Random.seed!(8101)
        local_rho = shadow(operations, 47; qubits=qubits)
        selected = qubits isa Integer ? [qubits] : qubits
        expected = partial_trace(full, fill(2, 3), setdiff(1:3, selected))
        @test size(local_rho) == (2^length(selected), 2^length(selected))
        @test local_rho ≈ expected atol=1e-12
        @test ishermitian(local_rho)
        @test tr(local_rho) ≈ 1 atol=1e-12
    end

    @test_throws ArgumentError shadow(operations, 0)
    @test_throws ArgumentError shadow(operations, -1)
    for qubits in ([0], [4], [1, 1], [1.0])
        @test_throws ArgumentError shadow(operations, 1; qubits=qubits)
    end
    for options in (Options(readout_noise=Noise1("bit_flip", 1.0)),
                    Options(measurement_mitigate=true))
        @test_throws ArgumentError shadow([Op("I", 1)], 1; options=options)
    end

    # Both Y eigenstates detect an incorrect rotation or outcome-sign convention.
    for (phase, sign) in (("S", 1), ("SD", -1))
        Random.seed!(8102)
        estimate = shadow([Op("H", 1), Op(phase, 1)], 4000; qubits=1)
        exact = (gate.I + sign * gate.Y) / 2
        @test norm(estimate - exact) < 0.10
    end

    # Tensoring averaged one-qubit estimates would lose these Bell correlations.
    bell = [Op("H", 1), Op("CX", 1, 2), Op("I", 3)]
    Random.seed!(8103)
    estimate = shadow(bell, 8000; qubits=[1, 2])
    bell_state = ComplexF64[1, 0, 0, 1] / sqrt(2)
    @test norm(estimate - bell_state * bell_state') < 0.13

    # Each shot must prepare a fresh noise/measurement trajectory. Reusing one
    # randomly collapsed state would give a pure state instead of I/2.
    mixed_preparations = (
        ([Op("I", 1)], Options(noise=NoiseModel("bit_flip", 0.5))),
        ([Op("I", 1), OpQC("bit_flip", 0.5, 1)], Options()),
        ([Op("H", 1), Op("MZ", 1)], Options()),
    )
    for (mixed, options) in mixed_preparations
        Random.seed!(8104)
        @test norm(shadow(mixed, 4000; qubits=1, options=options) - gate.I / 2) < 0.10
    end

    # A single snapshot has eigenvalues 2 and -1; it is not a physical state.
    snapshot = shadow([Op("I", 1)], 1)
    @test eigvals(Hermitian(Matrix(snapshot))) ≈ [-1, 2]
end

@testset "MPS classical shadows" begin
    # A complex Bell pair on nonadjacent sites tests Y signs, correlations and
    # subsystem ordering against an independent statevector partial trace.
    operations = [Op("H", 1), Op("S", 1), Op("CX", 1, 3)]
    exact = partial_trace(apply(zero_state(3), operations), [1, 3])
    Random.seed!(8201)
    estimate = shadow(operations, 1500; qubits=[1, 3], backend=:mps,
                      cutoff=0.0, maxdim=8)
    @test norm(estimate - exact) < 0.25

    Random.seed!(8202)
    full = shadow(operations, 31; backend=:tensor, cutoff=0.0)
    Random.seed!(8202)
    local_rho = shadow(operations, 31; qubits=[3, 1], backend="mps", cutoff=0.0)
    @test local_rho ≈ partial_trace(full, fill(2, 3), [2]) atol=1e-12
    @test ishermitian(local_rho)
    @test tr(local_rho) ≈ 1 atol=1e-12
    @test_throws ArgumentError shadow(operations, 1; backend=:density_matrix)

    # Fresh preparation is required for both Kraus noise and measurements.
    for (mixed, options) in (([Op("I", 1)], Options(noise=NoiseModel("bit_flip", 0.5))),
                             ([Op("H", 1), Op("MZ", 1)], Options()))
        Random.seed!(8203)
        estimate = shadow(mixed, 700; qubits=1, backend=:tensor, cutoff=0.0, options=options)
        @test norm(estimate - gate.I / 2) < 0.18
    end

    # A 72-site statevector cannot be constructed. The callback also verifies
    # that preparation and truncation options reach the MPS execution path.
    calls = Ref(0)
    check_mps = Op("check MPS", function (state; cutoff, maxdim)
        @test state isa it.MPS
        @test cutoff == 1e-12
        @test maxdim == 4
        calls[] += 1
        return state
    end)
    large = [Op("X", 1), Op("X", 72), check_mps]
    estimate = shadow(large, 3; qubits=[1, 72], backend=:mps,
                      cutoff=1e-12, maxdim=4)
    @test calls[] == 3
    @test size(estimate) == (4, 4)
    @test ishermitian(estimate)
    @test tr(estimate) ≈ 1 atol=1e-12
end
