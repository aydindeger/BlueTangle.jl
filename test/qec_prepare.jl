@testset "Physical Steane stabilizer-state preparation" begin
    code = steane_code()
    for basis in (:zero, :plus)
        target = basis == :zero ? zero_state(1) : plus_state(1)
        result = steane_prepare(basis)
        @test abs2(dot(target, code.decode(result.state))) ≈ 1
        @test get_N(result.state) == 7
        @test getproperty.(result.trace, :index) == collect(1:length(result.trace))
        @test count(loc -> loc.stage == :data_prep, result.trace) == 7
        @test 2 <= length(result.history) <= 4

        # Physical initialization faults are included in the same trace.
        injected = loc -> loc.stage == :data_prep && loc.op.qubit == 2 ? Op("X", 2) : Op[]
        faulty = steane_prepare(basis; fault=injected)
        diagnosed = qec_round(code, faulty.state; decoder=steane_decoder(code))
        @test abs2(dot(target, code.decode(diagnosed.state))) ≈ 1
    end
    mps = steane_prepare(:plus; backend=:tensor)
    @test mps.state isa it.MPS
    @test abs2(dot(plus_state(1), to_state(code.decode(mps.state)))) ≈ 1
    @test_throws ArgumentError steane_prepare(:magic)
    @test_throws ArgumentError steane_prepare(; backend=:invalid)
    @test_throws ArgumentError steane_prepare(; blocks=0)
    @test_throws ArgumentError steane_prepare(; prep_error=-0.1)
end
