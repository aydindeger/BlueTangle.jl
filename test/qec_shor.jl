@testset "Verified noisy Shor extraction" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.2+0.7im])
    encoded = code.encode(logical)
    schedule = shor_syndrome_circuit(code, 1)
    @test length(schedule) == 24
    @test count(e -> e.kind == :prep, schedule) == 5
    @test count(e -> e.kind == :measure, schedule) == 5
    @test [(e.op.qubit, e.op.target_qubit) for e in schedule if e.stage == :cat_verify && e.kind == :gate] == [(8, 12), (9, 12)]

    clean = shor_qec_round(code, encoded)
    @test get_N(clean.state) == code.n
    @test clean.history == [zeros(Int, code.m), zeros(Int, code.m)]
    @test length(clean.trace) == 288
    @test abs2(dot(logical, code.decode(clean.state))) ≈ 1
    @test all(event.index == i for (i, event) in enumerate(clean.trace))

    # The complete extraction acts on arbitrary complex logical amplitudes,
    # rather than testing only an encoded computational-basis state.
    for q in 1:code.n, name in ("X", "Y", "Z")
        damaged = apply(encoded, Op(name, q))
        result = shor_qec_round(code, damaged)
        @test result.syndrome == pauli_syndrome(code, Op(name, q))
        @test abs2(dot(logical, code.decode(result.state))) ≈ 1
        @test any(event -> event.stage == :recovery, result.trace)
    end

    # A dangerous two-cat X pattern from a *single* preparation-gate fault is
    # caught before interacting with data. Only the rejected ancilla is retried.
    hook = event -> event.round == 1 && event.check == 1 && event.attempt == 1 &&
                    event.stage == :cat_encode && event.op.name == "CX" &&
                    event.op.target_qubit == code.n+3 ? Op[Op("X", code.n+1)] : Op[]
    rejected = shor_qec_round(code, encoded; fault=hook)
    @test rejected.attempts[1][1] == 2
    @test length(rejected.trace) == 300
    @test abs2(dot(logical, code.decode(rejected.state))) ≈ 1
    @test_throws QECProtocolFailure shor_qec_round(code, encoded; fault=hook, max_attempts=1)

    # An erroneous check result must not immediately trigger data recovery.
    bad_record = event -> event.round == 1 && event.check == 1 &&
                          event.kind == :measure && event.stage == :cat_readout &&
                          event.op.qubit == code.n+1 ? Op[Op("X", code.n+1)] : Op[]
    repeated = shor_qec_round(code, encoded; fault=bad_record)
    @test length(repeated.history) == 3
    @test repeated.history[1] == [1, 0, 0, 0, 0, 0]
    @test all(iszero, repeated.history[end])
    @test abs2(dot(logical, code.decode(repeated.state))) ≈ 1
    @test_throws QECProtocolFailure shor_qec_round(code, encoded; fault=bad_record, max_rounds=2)

    # Recovery is a real noisy physical location, not a silent perfect frame
    # update. Injecting Z there leaves an error that a subsequent round removes.
    recovery_fault = event -> event.stage == :recovery ? Op[Op("Z", event.op.qubit)] : Op[]
    result = shor_qec_round(code, apply(encoded, Op("X", 1)); fault=recovery_fault)
    residual_syndrome = [round(Int, (1-expect(result.state, string_to_matrix(s)))/2) for s in code.stabilizers]
    @test residual_syndrome == pauli_syndrome(code, Op("Z", 1))
    repaired = shor_qec_round(code, result.state)
    @test abs2(dot(logical, code.decode(repaired.state))) ≈ 1

    # All zero-rate noise paths must reproduce the circuit without noise.
    result = shor_qec_round(code, encoded; noise=NoiseModel("depolarizing", 0.0),
                             prep_error=0.0, measurement_error=0.0)
    @test abs2(dot(logical, code.decode(result.state))) ≈ 1
    @test_throws QECProtocolFailure shor_qec_round(code, encoded; measurement_error=1.0, max_attempts=1)

    @test_throws ArgumentError shor_qec_round(code, encoded; prep_error=-0.1)
    @test_throws ArgumentError shor_qec_round(code, encoded; measurement_error=NaN)
    @test_throws ArgumentError shor_qec_round(code, encoded; max_rounds=1)
    @test_throws ArgumentError shor_qec_round(code, encoded; max_attempts=0)
    @test_throws DimensionMismatch shor_qec_round(code, zero_state(3))
    @test_throws ArgumentError shor_syndrome_circuit(code, 0)
    @test_throws ArgumentError shor_qec_round(code, encoded; fault=event -> Op("X", 1))
end

@testset "Native MPS Shor trajectory" begin
    code = steane_code()
    logical = normalize(ComplexF64[1, 0.3+0.4im])
    damaged = apply(code.encode(logical), Op("Y", 4))
    initial = to_MPS(damaged, N_MPS(code.n); cutoff=0.0)
    result = shor_qec_round(code, initial)
    @test result.state isa it.MPS
    @test get_N(result.state) == code.n
    @test result.syndrome == pauli_syndrome(code, Op("Y", 4))
    @test isapprox(abs2(inner(logical, code.decode(result.state))), 1; atol=1e-10)
    @test isapprox(abs2(inner(initial, damaged)), 1; atol=1e-10)
end
