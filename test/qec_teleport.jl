@testset "Steane logical teleportation" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.3+0.7im])
    encoded = code.encode(logical)
    branches = Dict{Tuple{Int,Int},Int}()
    audited_faults = 0

    # A complex superposition distinguishes a missing X or Z correction, even
    # though the logical Z probabilities alone are insensitive to a missing Z.
    for seed in 1:64
        Random.seed!(seed)
        result = steane_teleport(code, encoded)
        x, z = result.x_readout, result.z_readout
        @test abs2(dot(encoded, result.state)) ≈ 1
        @test x.corrected_qubit == z.corrected_qubit == 0
        @test count(e -> e.stage == :x_correction, result.trace) == 7 * z.outcome
        @test count(e -> e.stage == :z_correction, result.trace) == 7 * x.outcome
        @test [e.index for e in result.trace] == eachindex(result.trace)
        branch = (x.outcome, z.outcome)
        haskey(branches, branch) && continue
        branches[branch] = seed
        # Propagate every location Pauli through the sampled branch. After
        # its MZ, a measured wire has no more gates, so its final X component
        # records exactly whether its classical bit flipped.
        steps = filter(e -> e.component == :teleport, result.trace)
        paulis = ("I", "X", "Y", "Z")
        for (i, step) in enumerate(steps)
            wires = step.op.q == 1 ? (step.op.qubit,) : (step.op.qubit, step.op.target_qubit)
            later = Op[e.op for e in steps[i+1:end] if e.kind == :gate]
            for names in Iterators.product(ntuple(_ -> paulis, length(wires))...)
                all(==("I"), names) && continue
                errors = [Op(name, q) for (name, q) in zip(names, wires)]
                residual = propagate_pauli(errors, later; nqubits=21)
                flips = [e.qubit for e in residual if e.name in ("X", "Y")]
                @test count(q -> q in 1:7, flips) <= 1
                @test count(q -> q in 8:14, flips) <= 1
                @test count(e -> e.qubit in 15:21, residual) <= 1
                audited_faults += 1
            end
        end
        length(branches) == 4 && break
    end
    @test length(branches) == 4
    @test audited_faults == 1176

    # Physical position three in each block has wire number 3, 10 or 17.
    # One correlated CNOT fault can affect both classical words, while a
    # receiver fault may remain on the output until the next recovery.
    for (stage, q, errors, seed, xsite, zsite) in (
        (:bell_cnot,     10, [Op("Y", 10), Op("Y", 17)], 1, 3, 3),
        (:teleport_cnot,  3, [Op("Z", 3), Op("X", 10)],  1, 3, 3),
        (:x_readout,      3, [Op("X", 3)],               1, 3, 0),
        (:z_readout,     10, [Op("X", 10)],              1, 0, 3),
        (:x_correction,  17, [Op("Y", 17)], branches[(0, 1)], 0, 0),
        (:z_correction,  17, [Op("Y", 17)], branches[(1, 0)], 0, 0),
    )
        at_fault(e) = e.stage == stage && e.op.qubit == q
        Random.seed!(seed)
        result = steane_teleport(code, encoded; fault=e -> at_fault(e) ? errors : Op[])
        @test count(at_fault, result.trace) == 1
        @test result.x_readout.corrected_qubit == xsite
        @test result.z_readout.corrected_qubit == zsite
        recovered = qec_round(code, result.state; decoder=steane_decoder(code)).state
        @test abs2(dot(encoded, recovered)) ≈ 1
    end

    result = steane_teleport(code, apply(encoded, Op("Y", 4)))
    @test result.x_readout.corrected_qubit == result.z_readout.corrected_qubit == 4
    @test abs2(dot(encoded, result.state)) ≈ 1

    # Force one failed verification during preparation of the receiving block.
    reject(e) = e.component == :receiver_prepare && e.stage == :cat_verify && e.kind == :measure ?
        Op("X", e.op.qubit) : Op[]
    @test_throws QECProtocolFailure steane_teleport(code, encoded; fault=reject, max_attempts=1)
    @test_throws DimensionMismatch steane_teleport(code, zero_state(3))
    @test_throws ArgumentError steane_teleport(code, encoded; block=2)
    @test_throws ArgumentError steane_teleport(code, encoded; measurement_error=-0.1)

    # The first block is entangled with an asymmetric reference. The overlap
    # checks preservation of coherence and return of the receiver in slot one.
    two = steane_code(2)
    pair = normalize(ComplexF64[1, -0.2im, 0.4+0.3im, 0.7])
    input = two.encode(to_MPS(pair, N_MPS(2)))
    result = steane_teleport(two, input; block=1)
    @test result.state isa it.MPS
    @test isapprox(abs2(dot(pair, to_state(two.decode(result.state)))), 1; atol=1e-8)
    @test to_state(two.decode(input)) ≈ pair
end
