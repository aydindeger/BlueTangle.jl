@testset "Raw magic-state preparation" begin
    for dagger in (false, true)
        target = normalize(ComplexF64[1, cis((dagger ? -1 : 1)*pi/4)])
        raw = raw_magic_state(; dagger)
        @test abs2(dot(target, raw.state)) ≈ 1
        @test getproperty.(raw.trace, :stage) == [:raw_prep, :raw_h, :raw_t]
        @test all(!loc.phase_error_applied for loc in raw.trace if haskey(loc, :phase_error_applied))
        bad = raw_magic_state(; dagger, phase_error=1.0)
        @test isapprox(abs2(dot(target, bad.state)), 0; atol=1e-10)
        @test bad.trace[end].phase_error_applied
    end
    @test_throws ArgumentError raw_magic_state(; phase_error=-0.1)
    @test_throws ArgumentError raw_magic_state(; prep_error=NaN)
end

@testset "Physical 15-to-1 distillation circuit" begin
    supports = distillation15_supports()
    @test length(supports) == 15
    @test all(first(s) == 1 for s in supports)
    @test supports[1] == [1, 5]
    @test supports[15] == [1, 2, 3, 4, 5]
    target = normalize(ComplexF64[1, cis(pi/4)])
    observed = Set{Int}()
    for dagger in (false, true), seed in 1:3
        Random.seed!(seed)
        result = distill15(; dagger)
        expected = dagger ? conj(target) : target
        @test isapprox(abs2(dot(expected, result.state)), 1; atol=1e-9)
        @test result.attempts == 1
        @test result.raw_states == 15
        @test result.syndromes == [zeros(Int, 4)]
        @test getproperty.(result.trace, :index) == collect(1:length(result.trace))
        union!(observed, (loc.outcome for loc in result.trace if loc.stage == :distill_inject_readout))
    end
    @test observed == Set([0, 1])

    # Inject Z at raw preparation locations, not via a closed-form error sampler.
    # Every one- or two-input defect must be detected by the actual circuit.
    raw_fault(indices; first_attempt=false) = loc ->
        loc.stage == :raw_t && loc.check in indices && (!first_attempt || loc.attempt == 1) ?
        Op("Z", loc.op.qubit; noisy=false) : Op[]
    for j in 1:15
        @test_throws QECProtocolFailure distill15(; fault=raw_fault([j]), max_attempts=1)
    end
    for j in 1:14, k in j+1:15
        @test_throws QECProtocolFailure distill15(; fault=raw_fault([j, k]), max_attempts=1)
    end

    # The 35 undetected weight-three patterns yield the orthogonal magic state.
    accepted_triples = [(j, k, l) for j in 1:13 for k in j+1:14 for l in k+1:15 if xor(xor(j,k),l) == 0]
    @test length(accepted_triples) == 35
    for triple in accepted_triples
        result = distill15(; fault=raw_fault(triple), max_attempts=1)
        @test result.syndromes == [zeros(Int, 4)]
        @test isapprox(abs2(dot(gate.Z * target, result.state)), 1; atol=1e-9)
    end
    for triple in ((1,2,4), (2,4,8), (3,5,7))
        @test_throws QECProtocolFailure distill15(; fault=raw_fault(triple), max_attempts=1)
    end

    retried = distill15(; fault=raw_fault([1]; first_attempt=true), max_attempts=2)
    @test retried.attempts == 2
    @test retried.raw_states == 30
    @test retried.syndromes == [[0,0,0,1], [0,0,0,0]]
    @test isapprox(abs2(dot(target, retried.state)), 1; atol=1e-9)

    # This physical reference has unprotected Clifford gates. An output-gate
    # fault is accepted and spoils the resource, demonstrating the noise floor.
    output_fault = loc -> loc.stage == :distill_output ? Op("Z", 1; noisy=false) : Op[]
    faulty_output = distill15(; fault=output_fault)
    @test faulty_output.attempts == 1
    @test isapprox(abs2(dot(target, faulty_output.state)), 0; atol=1e-9)
    noisy = distill15(; noise=NoiseModel("depolarizing", 0.001),
                      prep_error=0.001, measurement_error=0.001, raw_error=0.01)
    @test norm(noisy.state) ≈ 1
    @test noisy.raw_states == 15noisy.attempts
    @test count(loc -> loc.stage == :raw_t, noisy.trace) == noisy.raw_states
    @test_throws ArgumentError distill15(; max_attempts=0)
    @test_throws ArgumentError distill15(; raw_error=Inf)
end

@testset "Exact ideal-Clifford distillation reference" begin
    # Independent enumeration of all 2^15 subsets of the 15 nonzero columns.
    # Acceptance is column XOR zero; odd accepted parity is the output error.
    counts = zeros(Int, 2, 16)
    for mask in 0:2^15-1
        check = 0
        for j in 1:15
            (mask & (1 << (j-1))) != 0 && (check = xor(check, j))
        end
        check == 0 || continue
        weight = count_ones(mask)
        counts[isodd(weight) ? 2 : 1, weight+1] += 1
    end
    @test sum(counts) == 2048
    @test counts[2, 4] == 35
    @test all(iszero, counts[:, 2:3])
    for p in (0.0, 0.0001, 0.01, 0.1, 0.5, 1.0)
        masses = [sum(counts[row, w+1]*p^w*(1-p)^(15-w) for w in 0:15) for row in 1:2]
        reference = distillation15_statistics(p)
        @test reference.acceptance ≈ sum(masses)
        @test reference.output_error ≈ masses[2]/sum(masses)
        @test reference.mean_raw_states ≈ 15/reference.acceptance
    end
    @test isapprox(distillation15_statistics(1e-5).output_error / 1e-15, 35; rtol=1e-3)
    @test_throws ArgumentError distillation15_statistics(-0.1)
end
