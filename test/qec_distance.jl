@testset "Code constructors and measured distances" begin
    for code in (shor_code(), steane_code(), reed_muller_code())
        @test code.k == 1
        result = code_distance(code)
        @test result.exact
        @test result.distance == result.lower_bound == 3
        @test length(result.witness) == 3
        @test all(iszero, pauli_syndrome(code, result.witness))
        action = logical_pauli(code, result.witness)
        @test any(!iszero, action.x) || any(!iszero, action.z)
        @test result.method == :css_support
    end
    shor = shor_code()
    @test (shor.n, shor.m) == (9, 8)
    @test (code_distance(shor).dx, code_distance(shor).dz) == (3, 3)
    # A weight-two stabilizer does not make this degenerate code distance two.
    @test logical_pauli(shor, [Op("Z", 1), Op("Z", 2)]) == (x=[0], z=[0])
    @test code_distance(StabilizerCode(shor.stabilizers; d=99)).distance == 3

    rm = reed_muller_code()
    @test (rm.n, rm.m) == (15, 14)
    result = code_distance(rm; max_weight=7)
    @test (result.distance, result.dx, result.dz) == (3, 7, 3)
    @test result.stop_reason == :found
    @test result.x_witness == collect(1:7)
    @test result.z_witness == [1, 2, 3]
    bounded = code_distance(rm; max_weight=4)
    @test bounded.exact && bounded.distance == 3
    @test bounded.dx === nothing && bounded.dx_lower_bound == 5
    @test bounded.dz == bounded.dz_lower_bound == 3
    interrupted = code_distance(rm; max_candidates=242)
    @test interrupted.truncated && interrupted.distance == 3 && interrupted.exact
    @test interrupted.dx === nothing && interrupted.dx_lower_bound == 3

    # Verify actual code words and the T/T† convention, independently of the
    # distance calculation and of any stored gate metadata.
    logical = normalize(ComplexF64[1, 0.2+0.7im])
    encoded = rm.encode(logical)
    @test rm.decode(encoded) ≈ logical
    @test all(s -> isapprox(expect(encoded, string_to_matrix(s)), 1; atol=1e-10), rm.stabilizers)
    for name in ("T", "TD")
        transformed = rm.decode(rm.apply(encoded, name, 1))
        @test isapprox(transformed, gates(name) * logical; atol=1e-8)
    end
    all_x = apply(encoded, [Op("X", q) for q in 1:15])
    all_z = apply(encoded, [Op("Z", q) for q in 1:15])
    @test rm.decode(all_x) ≈ gate.X * logical
    @test rm.decode(all_z) ≈ gate.Z * logical
end

@testset "Bounded and non-CSS distance search" begin
    code = steane_code()
    result = code_distance(code; max_weight=2)
    @test !result.exact && result.distance === nothing && result.witness === nothing
    @test result.lower_bound == 3
    @test result.searched_through == 2
    @test result.stop_reason == :max_weight
    @test !result.truncated
    partial = code_distance(code; max_candidates=1)
    @test partial.truncated && partial.candidates == 1
    @test partial.lower_bound == 1 && partial.searched_through == 0
    frontier = code_distance(code; max_candidates=14)
    @test frontier.truncated && frontier.candidates == 14
    @test frontier.lower_bound == 2 && frontier.searched_through == 1

    five = StabilizerCode(["X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"]; d=3)
    general = code_distance(five)
    @test general.exact && general.distance == 3
    @test general.method == :pauli_support
    @test general.dx === general.dz === nothing
    @test general.dx_lower_bound === general.dz_lower_bound === nothing
    @test all(iszero, pauli_syndrome(five, general.witness))
    @test any(!iszero, vcat(logical_pauli(five, general.witness)...))
    partial_general = code_distance(five; max_candidates=15)
    @test partial_general.lower_bound == 2 && partial_general.searched_through == 1

    repetition = StabilizerCode(["Z,Z,I", "I,Z,Z"])
    @test code_distance(repetition).distance == 1 # phase errors are unprotected
    @test_throws ArgumentError code_distance(code; max_weight=0)
    @test_throws ArgumentError code_distance(code; max_candidates=0)
    @test_throws ArgumentError code_distance(StabilizerCode(["Z"]))
end
