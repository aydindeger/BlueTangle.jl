import PhasePoly

@testset "Optional PhasePoly backend" begin
    programs = [
        [Op("X",1),Op("I",3)],
        [Op("H",1),Op("H",3),Op("CX",1,2),Op("T",3),Op("Y",2),Op("CZ",2,3),Op("H",3)],
        [Op("H",1),Op("H",2),Op("P(pi/8)",1),Op("CP(pi/4)",1,2),Op("RZ(pi/4)",2;control=1),Op("I",3)],
        [Op("RZ(2pi)",1),Op("Y",3)],
        [Op("H",1),Op("T",1),Op("H",1),Op("X",3)],
        [Op("X",1),Op("H",3),Op("H",3)],
        [Op("H",1),Op("H",2),Op("H",3),Op("CCZ",1,2,3),Op("CTD",1,3)],
        [Op("H",1),Op("Y",2;control=1),Op("S",3;control=2)],
    ]
    for ops in programs
        expected = apply(zero_state(3),ops)
        compact = phasepoly(ops; nqubits=3)
        @test get_N(compact) == 3
        @test to_state(compact) ≈ expected atol=1e-9
        @test to_state(phasepoly(ops)) ≈ expected atol=1e-9
        for index in 0:7
            @test amplitude(compact,index) ≈ expected[index+1] atol=1e-9
            @test amplitude(compact,int2bin(index,3)) ≈ expected[index+1] atol=1e-9
        end
        for axis in ("X","Y","Z")
            @test expect(compact,axis) ≈ expect(expected,axis) atol=1e-9
            @test expect(compact,axis,1) ≈ expect(expected,Op(axis,1)) atol=1e-9
        end
        @test correlation(compact,"X,Z",[1,3]) ≈ correlation(expected,"X,Z",[1,3]) atol=1e-9
        @test expect(compact,"XIY") ≈ correlation(expected,"X,Y",[1,3]) atol=1e-9
    end

    # This query must use PhasePoly's compact path: 2^80 amplitudes cannot be materialized.
    large = phasepoly([Op("H",q) for q in 1:80])
    @test get_N(large) == 80
    @test expect(large,"X",71) ≈ 1
    @test correlation(large,"X,X",[2,79]) ≈ 1
    @test_throws ArgumentError phasepoly([Op("MZ",1)])
    @test_throws ArgumentError phasepoly([Op("RES",1)])
    @test_throws ArgumentError phasepoly([Op("RX(0.1)",1)])
    @test_throws ArgumentError phasepoly([Op("P(pi/16)",1)])
    @test_throws ArgumentError phasepoly([Op("X",gate.Z,1)])
    @test_throws ArgumentError phasepoly([Op("CCX",1,2,3)])

    # Unified backend=:phasepoly interface (apply, measure, sample, shadow, HSP/HY)
    ops_pp = [Op("H", 1), Op("H", 2), Op("CZ", 1, 2), Op("T", 1), Op("HSP", 2), Op("HY", 1)]
    sv_pp = apply(zero_state(2), ops_pp)
    st_pp = apply(ops_pp; backend=:phasepoly)
    @test to_state(st_pp) ≈ sv_pp atol=1e-9
    @test to_rho(st_pp) ≈ sv_pp * sv_pp' atol=1e-9
    @test fidelity(st_pp, sv_pp) ≈ 1.0 atol=1e-9
    @test measure(st_pp).expect ≈ expect(sv_pp, "Z") atol=1e-9
    @test to_state(phasepoly(ops_pp; N=2)) ≈ sv_pp atol=1e-9

    m_pp = measure(ops_pp, 2000; backend=:phasepoly)
    @test m_pp.expect ≈ expect(sv_pp, "Z") atol=0.10
    @test isempty(m_pp.mid_measurements)
    @test length(sample(ops_pp, 25; backend=:phasepoly)) == 25
    @test measure([Op("H", 1)], 20; backend=:phasepoly, options=Options(measurement_basis="X")).bitstr == [0]
    @test measure([Op("H", 1), Op("S", 1)], 20; backend=:phasepoly, options=Options(measurement_basis="Y")).bitstr == [0]

    ghz_pp = [Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
    rho_13 = shadow(ghz_pp, 2000; qubits=[1, 3], backend=:phasepoly)
    @test size(rho_13) == (4, 4)
    @test norm(Matrix(rho_13) - Matrix(partial_trace(apply(zero_state(3), ghz_pp), [1, 3]))) < 0.12

    @test_throws ArgumentError apply(ops_pp; backend=:phasepoly, noise=NoiseModel("depolarizing", 0.1))
    @test_throws ArgumentError measure(ops_pp, 10; backend=:phasepoly, options=Options(noise=NoiseModel("depolarizing", 0.1)))
    @test_throws ArgumentError measure(ops_pp, 10; backend=:phasepoly, options=Options(readout_noise=Noise1("bit_flip", 0.1)))
    @test_throws ArgumentError measure([Op("MZ", 1)], 1; backend=:phasepoly)

    # Tuple shorthand and mixed lists with phasepoly and backend=:phasepoly
    tup_pp = [("H", 1), ("CX", 1, 2), ("T", 1)]
    mix_pp = [("H", 1), Op("CX", 1, 2), ("T", 1)]
    sv_tup = apply(tup_pp)
    @test to_state(phasepoly(tup_pp)) ≈ sv_tup atol=1e-9
    @test to_state(phasepoly(mix_pp)) ≈ sv_tup atol=1e-9
    @test to_state(apply(tup_pp; backend=:phasepoly)) ≈ sv_tup atol=1e-9
    @test to_state(apply(mix_pp; backend=:phasepoly)) ≈ sv_tup atol=1e-9
    @test all(b -> b in (0, 3), measure(tup_pp, 50; backend=:phasepoly).bitstr)
    @test all(b -> b in (0, 3), sample(tup_pp, 25; backend=:phasepoly))
    @test size(shadow(tup_pp, 20; qubits=1, backend=:phasepoly)) == (2, 2)
end
