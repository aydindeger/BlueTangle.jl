@testset "Composable ideal QEC rounds" begin
    code = StabilizerCode([
        "X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"
    ]; d=3)
    logical = normalize(ComplexF64[1, 0.2+0.7im])
    encoded = code.encode(logical)

    # A round leaves only data wires, so several independent error/recovery
    # cycles compose without decoding or accumulating syndrome ancillas.
    state = copy(encoded)
    for op in (Op("Y", 2), Op("X", 5), Op("Z", 3))
        state = apply(state, op)
        result = qec_round(code, state)
        expected = [round(Int, (1-expect(state, string_to_matrix(s)))/2) for s in code.stabilizers]
        @test result.syndrome == expected
        @test get_N(result.state) == code.n
        @test abs2(dot(logical, code.decode(result.state))) ≈ 1
        state = result.state
    end
    @test qec_round(code, state).syndrome == zeros(Int, code.m)

    # Projective syndrome readout resolves a coherent superposition of two
    # errors. Either record must recover the same logical amplitudes.
    damaged = normalize(apply(encoded, Op("X", 1)) + im * apply(encoded, Op("Z", 3)))
    mps = to_MPS(damaged, N_MPS(code.n); cutoff=0.0)
    for seed in (1401, 1402)
        Random.seed!(seed)
        dense_result = qec_round(code, damaged)
        Random.seed!(seed)
        mps_result = qec_round(code, mps)
        @test mps_result.state isa it.MPS
        @test mps_result.syndrome == dense_result.syndrome
        @test abs2(inner(mps_result.state, dense_result.state)) ≈ 1
        @test abs2(dot(logical, code.decode(dense_result.state))) ≈ 1
        @test get_N(qec_round(code, mps_result.state).state) == code.n
    end
    @test to_state(mps) ≈ damaged # input MPS was not mutated

    logical_rho = ComplexF64[0.6 0.2+0.1im; 0.2-0.1im 0.4]
    rho = code.encode(logical_rho)
    noisy_rho = 0.4rho + 0.3apply(rho, Op("X", 1)) + 0.3apply(rho, Op("Y", 4))
    result = qec_round(code, noisy_rho)
    @test result.syndrome === nothing
    @test size(result.state) == size(rho)
    @test tr(result.state) ≈ tr(noisy_rho)
    @test code.decode(result.state) ≈ logical_rho
    @test qec_round(code, result.state).state ≈ result.state
    @test code.decode(qec_round(code, damaged * damaged').state) ≈ logical * logical'

    # The same callback contract permits an error-model-specific decoder.
    repetition = StabilizerCode(["Z,Z,I,I,I", "I,Z,Z,I,I", "I,I,Z,Z,I", "I,I,I,Z,Z"])
    recovery = [Op("X", 1; noisy=false), Op("X", 2; noisy=false)]
    decoder = bits -> all(iszero, bits) ? Op[] : recovery
    encoded_rep = repetition.encode(logical)
    result = qec_round(repetition, apply(encoded_rep, recovery); decoder=decoder)
    @test result.syndrome == [0, 1, 0, 0]
    @test abs2(dot(logical, repetition.decode(result.state))) ≈ 1
    @test_throws ArgumentError qec_round(repetition, apply(encoded_rep, recovery); decoder=_ -> Op[])

    @test_throws DimensionMismatch qec_round(code, zero_state(3))
    @test_throws DimensionMismatch qec_round(code, to_rho(zero_state(3)))
    @test_throws ArgumentError qec_round(code, zeros(ComplexF64, 2^code.n))
    @test_throws ArgumentError qec_round(code, zeros(ComplexF64, 2^code.n, 2^code.n))
end
