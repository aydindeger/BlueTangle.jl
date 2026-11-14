@testset "Sampling and measurement statistics" begin
    amplitudes = ComplexF64[0, 1 / sqrt(3), 0, im * sqrt(2 / 3)]
    expected = ([1, 3], [1 / 3, 2 / 3])
    for state in (amplitudes, view(amplitudes, :), sa.sparse(amplitudes))
        outcomes, probabilities = sample_exact(state)
        @test outcomes == expected[1]
        @test probabilities ≈ expected[2]
        @test measure(state).expect ≈ [1 / 3 - 2 / 3, -1]
    end

    # A sparse 20-qubit state should not need a full Hilbert-space weight table.
    last_state = sa.sparsevec([2^20], [1.0], 2^20)
    @test sample(last_state, 100) == fill(2^20 - 1, 100)
    @test isempty(sample(last_state, 0))
    draws = sample(sa.sparse(amplitudes), 10_000)
    @test all(in((1, 3)), draws)
    @test isapprox(count(==(1), draws) / length(draws), 1 / 3; atol=0.03)

    probabilities = fill(1 / 2^6, 2^6)
    outcomes = 0:2^6-1
    @test mag_moments(6, outcomes, probabilities, 1) == 0
    @test mag_moments(6, outcomes, probabilities, 2) == 6
    @test mag_moments(6, outcomes, probabilities, 4) == 3 * 6^2 - 2 * 6
    @test mag_moments(40, [0], [1.0], 12) == 40.0^12
    @test mag_moments(40, [2^40-1], [1.0], 11) == (-40.0)^11
    @test mag_moments(6, collect(outcomes), sa.sparse(probabilities), 2) == 6
    @test mag_moments(0, [0], [1.0], 2) == 0

    focks = int2bin.(outcomes, 6)
    @test BlueTangle._sample_to_expectation((focks, probabilities), [1, 6]) == 0
    @test BlueTangle._sample_to_expectation((focks, probabilities), [2, 2]) == 1
    @test BlueTangle._sample_to_expectation((focks, probabilities), Int[]) == 1
end
