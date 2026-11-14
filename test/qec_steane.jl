@testset "Steane logical Clifford circuits" begin
    @test_throws ArgumentError steane_code(0)
    code = steane_code(2)
    @test (code.n, code.k, code.d, code.m) == (14, 2, 3, 12)
    @test code.logical_indices == [7, 14]
    # A complex entangled input detects wrong S/SD conventions and block order.
    logical = normalize(ComplexF64[1, 0.3+0.2im, -0.6im, -0.4+0.1im])
    encoded = code.encode(logical)
    @test code.decode(encoded) ≈ logical
    for q in 1:2, name in ("I", "X", "Y", "Z", "H", "S", "SD")
        actual = code.decode(code.apply(encoded, name, q))
        @test actual ≈ apply(logical, Op(name, q))
    end
    for (q, t) in ((1, 2), (2, 1)), name in ("CX", "CNOT", "CZ", "SWAP")
        actual = code.decode(code.apply(encoded, name, q, t))
        @test actual ≈ apply(logical, Op(name, q, t))
        @test all(op -> div(op.qubit-1, 7) != div(op.target_qubit-1, 7), code.ops(name, q, t))
    end
    algorithm = [Op("H", 2), Op("S", 1), Op("CX", 2, 1), Op("SD", 2)]
    @test code.decode(apply(code.ops(algorithm), encoded)) ≈ apply(algorithm, logical)
    @test_throws KeyError code.ops("T", 1)
    # Logical compilation must not silently drop a control or a custom matrix.
    for invalid in (Op("X", 1; control=2), Op("H", gate.S, 1), Op("X", gate.Z, 1))
        @test_throws ArgumentError code.ops(invalid)
        @test_throws ArgumentError code.ops([invalid])
        @test_throws ArgumentError code.apply(encoded, invalid)
    end
    custom = StabilizerCode(BlueTangle._STEANE_CHECKS, Dict(("custom_H", 1) => code.ops("H", 1)))
    @test custom.decode(custom.apply(custom.encode(), Op("custom_H", gate.H, 1))) ≈ plus_state(1)

    # Density and MPS paths exercise the same logical compilation on one block.
    single = steane_code()
    psi = normalize(ComplexF64[1, 0.2+0.5im])
    gates_logical = [Op("H", 1), Op("S", 1), Op("Y", 1)]
    target = apply(gates_logical, psi)
    rho = single.encode(psi * psi')
    @test single.decode(apply(single.ops(gates_logical), rho)) ≈ target * target'
    mps = single.encode(to_MPS(psi, N_MPS(1); cutoff=0.0))
    @test to_state(single.decode(apply(single.ops(gates_logical), mps))) ≈ target
end

@testset "Steane block decoder" begin
    code = steane_code(2)
    decoder = steane_decoder(code)
    @test isempty(decoder(zeros(Int, code.m)))
    @test_throws DimensionMismatch decoder([0])
    @test_throws ArgumentError decoder(fill(2, code.m))
    @test_throws ArgumentError steane_decoder(StabilizerCode(["Z,Z,I", "I,Z,Z"]))
    reordered = StabilizerCode(reverse(BlueTangle._STEANE_CHECKS))
    @test_throws ArgumentError steane_decoder(reordered)

    # Independently calculate check anticommutation; every combination of one
    # Pauli per block must decode to precisely those physical errors.
    for q in 1:7, t in 8:14, a in ("X", "Y", "Z"), b in ("X", "Y", "Z")
        errors = [Op(a, q), Op(b, t)]
        bits = zeros(Int, code.m)
        for (i, row) in enumerate(code.stabilizers)
            factors = split(row, ",")
            bits[i] = mod(sum(factors[e.qubit] != "I" && factors[e.qubit] != e.name for e in errors), 2)
        end
        recovery = decoder(bits)
        @test [(e.name, e.qubit) for e in recovery] == [(a, q), (b, t)]
    end
end
