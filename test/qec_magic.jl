@testset "Steane destructive logical readout" begin
    checks = stabilizers_to_generator(BlueTangle._STEANE_CHECKS)[1:3, 1:7]
    # Exhaust all classical Hamming codewords and all single readout flips.
    words = [int2bin(i, 7) for i in 0:127 if all(iszero, mod.(checks * int2bin(i, 7), 2))]
    @test length(words) == 16
    for word in words, q in 0:7
        corrupted = copy(word)
        q > 0 && (corrupted[q] ⊻= 1)
        result = steane_readout(corrupted)
        @test result.corrected_bits == word
        @test result.outcome == mod(sum(word), 2)
        @test result.corrected_qubit == q
    end
    @test_throws DimensionMismatch steane_readout([0, 1])
    @test_throws ArgumentError steane_readout(fill(2, 7))

    code = steane_code()
    for (basis, logical, expected) in ((:Z, zero_state(1), 0), (:Z, one_state(1), 1),
                                      (:X, plus_state(1), 0), (:X, minus_state(1), 1))
        encoded = code.encode(logical)
        for q in 1:7, name in ("X", "Y", "Z")
            result = steane_measure(code, apply(encoded, Op(name, q)); basis=basis)
            @test result.outcome == expected
            @test length(result.state) == 1
            @test norm(result.state) ≈ 1
        end
    end
    # Destructively measuring one half of an encoded Bell pair keeps the
    # remaining block encoded and conditionally correlated with the result.
    two = steane_code(2)
    bell = apply(zero_state(2), [Op("H", 1), Op("CX", 1, 2)])
    result = steane_measure(two, two.encode(bell); block=1)
    @test get_N(result.state) == 7
    @test abs2(dot(product_state([result.outcome]), code.decode(result.state))) ≈ 1
    @test_throws ArgumentError steane_measure(code, code.encode(); basis=:Y)
    @test_throws ArgumentError steane_measure(code, code.encode(); block=2)
end

@testset "Encoded T and T-dagger injection" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.3+0.7im])
    encoded = code.encode(logical)
    phase = ComplexF64[1 0; 0 cis(pi/4)]
    resources = Dict(false => code.encode(phase * plus_state(1)),
                     true => code.encode(phase' * plus_state(1)))
    observed = Set{Int}()
    for dagger in (false, true), seed in 1:8
        Random.seed!(seed)
        result = steane_t(code, encoded, resources[dagger]; dagger=dagger)
        expected = (dagger ? phase' : phase) * logical
        @test isapprox(abs2(dot(expected, code.decode(result.state))), 1; atol=1e-10)
        @test get_N(result.state) == 7
        @test length(result.trace) == (result.outcome == 0 ? 14 : 21)
        push!(observed, result.outcome)
    end
    @test observed == Set([0, 1])

    # Applying T to either block preserves arbitrary entanglement with the other.
    two = steane_code(2)
    pair = normalize(ComplexF64[1, 0.2im, -0.4+0.3im, 0.7])
    for block in 1:2
        result = steane_t(two, two.encode(pair), resources[false]; block=block)
        expected = apply(pair, Op("exact T", phase, block))
        @test isapprox(abs2(dot(expected, two.decode(result.state))), 1; atol=1e-10)
    end

    # Resource logical faults survive stabilizer correction: magic quality is
    # a required input, not something the consumption gadget can verify.
    bad_magic = code.apply(resources[false], "Z", 1)
    result = steane_t(code, code.encode(:plus), bad_magic)
    recovered = code.decode(qec_round(code, result.state).state)
    target = phase * plus_state(1)
    @test isapprox(abs2(dot(target, recovered)), 0; atol=1e-10)
    @test isapprox(abs2(dot(gate.Z * target, recovered)), 1; atol=1e-10)

    @test_throws DimensionMismatch steane_t(code, zero_state(3), resources[false])
    @test_throws DimensionMismatch steane_t(code, encoded, zero_state(1))
    @test_throws ArgumentError steane_t(code, encoded, resources[false]; measurement_error=-0.1)
end

@testset "Single-fault containment in encoded T consumption" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.3+0.7im])
    phase = ComplexF64[1 0; 0 cis(pi/4)]
    encoded = code.encode(logical)
    magic = code.encode(phase * plus_state(1))
    expected = phase * logical
    paulis = ("I", "X", "Y", "Z")
    # One run in each logical measurement branch. The seeds are discovered
    # from the actual sampler rather than depending on Julia RNG constants.
    branch_seeds = Dict{Int,Int}()
    for seed in 1:32
        Random.seed!(seed)
        outcome = steane_t(code, encoded, magic).outcome
        get!(branch_seeds, outcome, seed)
        length(branch_seeds) == 2 && break
    end
    @test length(branch_seeds) == 2
    fidelity_after_ideal_ec(state) = abs2(dot(expected, code.decode(qec_round(code, state).state)))
    fault_outcomes = Dict{Tuple{Int,String,String},Set{Int}}()
    for seed in values(branch_seeds), q in 1:7, a in paulis, b in paulis
        a == "I" && b == "I" && continue
        fault = loc -> loc.stage == :t_cnot && loc.op.qubit == q ?
            [Op(a, q; noisy=false), Op(b, code.n+q; noisy=false)] : Op[]
        Random.seed!(seed)
        result = steane_t(code, encoded, magic; fault=fault)
        push!(get!(fault_outcomes, (q, a, b), Set{Int}()), result.outcome)
        @test isapprox(fidelity_after_ideal_ec(result.state), 1; atol=1e-10)
    end
    @test length(fault_outcomes) == 105
    @test all(outcomes == Set([0, 1]) for outcomes in values(fault_outcomes))
    for seed in values(branch_seeds), q in 1:7, name in ("X", "Y", "Z")
        Random.seed!(seed)
        result = steane_t(code, apply(encoded, Op(name, q)), magic)
        @test isapprox(fidelity_after_ideal_ec(result.state), 1; atol=1e-10)
        Random.seed!(seed)
        result = steane_t(code, encoded, apply(magic, Op(name, q)))
        @test isapprox(fidelity_after_ideal_ec(result.state), 1; atol=1e-10)
    end
    for seed in values(branch_seeds), q in 1:7
        fault = loc -> loc.stage == :t_readout && loc.op.qubit == code.n+q ?
            Op("X", code.n+q; noisy=false) : Op[]
        Random.seed!(seed)
        result = steane_t(code, encoded, magic; fault=fault)
        @test isapprox(fidelity_after_ideal_ec(result.state), 1; atol=1e-10)
    end
    for q in 1:7, name in ("X", "Y", "Z")
        fault = loc -> loc.stage == :t_correction && loc.op.qubit == q ?
            Op(name, q; noisy=false) : Op[]
        Random.seed!(branch_seeds[1])
        result = steane_t(code, encoded, magic; fault=fault)
        @test result.outcome == 1
        @test isapprox(fidelity_after_ideal_ec(result.state), 1; atol=1e-10)
    end
    Random.seed!(149)
    result = steane_t(code, encoded, magic; noise=NoiseModel("depolarizing", 0.001),
                      measurement_error=0.001)
    @test norm(result.state) ≈ 1
    @test get_N(result.state) == 7
end

@testset "Native MPS magic consumption" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.3+0.7im])
    phase = ComplexF64[1 0; 0 cis(pi/4)]
    # Deliberately reuse the same site Index identities between the two blocks.
    sites = N_MPS(7)
    data = to_MPS(code.encode(logical), sites; cutoff=0.0)
    magic = to_MPS(code.encode(phase * plus_state(1)), sites; cutoff=0.0)
    saved_data, saved_magic = to_state(data), to_state(magic)
    result = steane_t(code, data, magic)
    @test result.state isa it.MPS
    @test isapprox(abs2(dot(phase*logical, to_state(code.decode(result.state)))), 1; atol=1e-8)
    @test to_state(data) ≈ saved_data
    @test to_state(magic) ≈ saved_magic
    measurement = steane_measure(code, result.state; basis=:X)
    @test length(measurement.state) == 1
    @test norm(measurement.state) ≈ 1
    @test_throws ArgumentError steane_t(code, data, saved_magic)
end
