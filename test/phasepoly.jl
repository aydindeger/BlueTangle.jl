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
end
