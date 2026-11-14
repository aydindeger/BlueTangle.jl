@testset "Raw encoded magic resources" begin
    code = steane_code()
    for dagger in (false, true)
        target = normalize(ComplexF64[1, cis((dagger ? -1 : 1)*pi/4)])
        result = steane_raw_magic(; dagger)
        @test abs2(dot(target, code.decode(result.state))) ≈ 1
        @test !result.raw_phase_error
        @test getproperty.(result.trace, :index) == collect(1:length(result.trace))
        bad = steane_raw_magic(; dagger, raw_error=1.0)
        @test abs2(dot(gate.Z * target, code.decode(bad.state))) ≈ 1
        @test bad.raw_phase_error
    end
    @test_throws ArgumentError steane_raw_magic(; raw_error=NaN)
    @test_throws ArgumentError steane_distill15(; max_attempts=0)
    @test_throws ArgumentError steane_distill15(; cutoff=-1)
    # Gate definitions must retain precision through fifteen rotations.
    @test norm(gate.T' * gate.T - I) < 1e-15
    @test norm(gate.TD' * gate.TD - I) < 1e-15
    @test norm(gate.CT' * gate.CT - I) < 1e-15
    @test norm(gate.CTD' * gate.CTD - I) < 1e-15
end

@testset "Removing measured MPS wires" begin
    state = apply(zero_state(5), [Op("H", 2), Op("CX", 2, 4), Op("S", 4)])
    mps = to_MPS(state, N_MPS(5); cutoff=0.0)
    for keep in ([2, 4], [2, 3, 4, 5], [1, 2, 3, 4], [4, 2])
        @test to_state(reduce_and_relabel_qubits(mps, keep)) ≈ reduce_and_relabel_qubits(state, keep)
    end
    @test_throws ArgumentError reduce_and_relabel_qubits(mps, [2])
    @test only(reduce_and_relabel_qubits(zero_state(N_MPS(3)), Int[])) ≈ 1
end

@testset "Encoded physical 15-to-1 factory" begin
    code = steane_code()
    target = normalize(ComplexF64[1, cis(pi/4)])
    first = steane_distill15(; max_attempts=1, cutoff=1e-12)
    @test abs2(dot(target, code.decode(first.state))) ≈ 1
    @test first.raw_states == 15
    @test first.syndromes == [zeros(Int, 4)]
    @test getproperty.(first.trace, :index) == collect(1:length(first.trace))
    @test count(e -> e.component == :raw && e.stage == :raw_t, first.trace) == 15
    @test count(e -> e.component == :parity, first.trace) == 16*7
    @test first.history[end].accepted

    # A fault at the output phase layer is repaired by the actual final Shor
    # circuit. Tightening cutoff also checks the numerical result.
    fault = loc -> loc.component == :output_phase && loc.op.qubit == 1 ? Op("Z", 1) : Op[]
    second = steane_distill15(; max_attempts=1, cutoff=1e-14, fault=fault)
    @test abs2(dot(target, code.decode(second.state))) ≈ 1
    @test abs2(dot(first.state, second.state)) ≈ 1
    @test any(e -> e.component == :output_correction && e.stage == :recovery, second.trace)
    @test get_N(second.state) == 7
    consumed = steane_t(code, steane_prepare(:plus).state, first.state)
    actual = shor_qec_round(code, consumed.state)
    @test abs2(dot(target, code.decode(actual.state))) ≈ 1

    # One bad raw input is rejected; all its resources still count. The next
    # attempt creates fresh states instead of replaying a sampled preparation.
    raw_fault = loc -> loc.factory_attempt == 1 && loc.component == :raw &&
        loc.resource == 1 && loc.stage == :raw_t ? Op("Z", loc.op.qubit) : Op[]
    retried = steane_distill15(; dagger=true, max_attempts=2, fault=raw_fault)
    @test retried.attempts == 2
    @test retried.raw_states == 30
    @test getproperty.(retried.history, :accepted) == [false, true]
    @test abs2(dot(conj(target), code.decode(retried.state))) ≈ 1
end
