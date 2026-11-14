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

@testset "Clear measurement and execution results" begin
    # Integer amplitudes are a quantum state, not a list of previously drawn shots.
    excited = measure([0, 1], 8)
    @test excited.number_of_qubits == 1
    @test excited.number_of_experiment == 8
    @test excited.bitstr == [1]
    @test excited.sample == [1.0]
    histogram = Measurement([0, 1, 1, 1], 1; label="recorded shots")
    @test histogram.sample == [0.25, 0.75]
    @test histogram.expect == [-0.5]
    @test histogram.circuit_name == "recorded shots"

    operations = [Op("H", 1), Op("CX", 1, 2)]
    for backend in (:statevector, :density_matrix, :mps)
        Random.seed!(317)
        draws = sample(operations, 24; backend=backend)
        @test length(draws) == 24
        @test all(in((0, 3)), draws)
        Random.seed!(317)
        summary = measure(operations, 24; backend=backend)
        outcomes, frequencies = get_probs_from_sample(draws, 2)
        @test summary.bitstr == outcomes
        @test summary.sample == frequencies

        # A record describes explicit measurements, not the final state after reset.
        reset = [Op("X", 1), Op("MZ", 1), Op("RESET", 1)]
        @test run(reset, 3; backend=backend) == [[1], [1], [1]]
        @test sample(reset, 3; backend=backend) == [0, 0, 0]
        @test sample(Op[], 2; N=2, backend=backend) == [0, 0]
        @test isempty(sample(operations, 0; backend=backend))
        @test sample([Op("I", 1)], 3; backend=backend,
                     options=Options(readout_noise=Noise1("bit_flip", 1.0))) == [1, 1, 1]
    end
    @test_throws ArgumentError sample(operations, 3; options=Options(measurement_mitigate=true))
    @test_throws ArgumentError sample(operations, 3; options=Options(density_matrix=true))
    @test_throws ArgumentError sample(operations, -1)
    @test_throws ArgumentError measure(operations, 0)

    plus = plus_state(1)
    rho = to_rho(plus)
    @test apply(rho, Op("MZ", 1)) ≈ [0.5 0; 0 0.5]
    for initial in (plus, rho, to_MPS(plus, N_MPS(1)))
        before = initial isa it.MPS ? to_state(initial) : copy(initial)
        result = run(initial, Op("MZ", 1))
        @test keys(result) == (:state, :record)
        @test result.record in ([0], [1])
        final = result.state isa it.MPS ? to_state(result.state) : result.state
        final_rho = final isa AbstractMatrix ? final : to_rho(final)
        expected = only(result.record) == 0 ? [1 0; 0 0] : [0 0; 0 1]
        @test final_rho ≈ expected
        @test (initial isa it.MPS ? to_state(initial) : initial) ≈ before
        @test measure(initial).sample ≈ [0.5, 0.5]
    end
    result = run(one_state(1), [("MZ", 1), ("RESET", 1)])
    @test result.record == [1]
    @test result.state ≈ zero_state(1)
end
