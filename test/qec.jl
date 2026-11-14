@testset "Stabilizer encoding and correction" begin
    # Check the phase rules independently against dense Pauli conjugation.
    pauli_names = ("I", "X", "Z", "Y")
    for name in ("H", "S", "X", "Y", "Z", "CX", "CY", "CZ")
        op = name in ("CX", "CY", "CZ") ? Op(name, 1, 2) : Op(name, 1)
        U = Matrix(op.expand(2))
        for p1 in 1:4, p2 in 1:4
            v = [mod(p1-1, 2), mod(p2-1, 2), (p1-1) ÷ 2, (p2-1) ÷ 2]
            P = kron(gates(pauli_names[p1]), gates(pauli_names[p2]))
            for inverse in (false, true)
                result, sign = BlueTangle._qec_conjugate_pauli(v, [op]; inverse=inverse)
                names = [pauli_names[1+result[q]+2result[2+q]] for q in 1:2]
                expected = inverse ? U' * P * U : U * P * U'
                @test sign * kron(gates.(names)...) ≈ expected
            end
        end
    end
    @test BlueTangle.logical_vec_2_ops([1, 1]) == [("Y", 1)]

    code_rows = [
        ["Z,Z,I", "I,Z,Z"],
        ["X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"],
        ["X,X,X,X,I,I,I", "X,X,I,I,X,X,I", "X,I,X,I,X,I,X",
         "Z,Z,Z,Z,I,I,I", "Z,Z,I,I,Z,Z,I", "Z,I,Z,I,Z,I,Z"],
        ["Y,Y"],
        ["I,X,X", "Z,Z,Z"],
        # Binary row operations create a negative Pauli phase in this code.
        ["Y,X,Y", "Y,Z,X"],
        ["X,X,X,X", "Z,Z,Z,Z"],
    ]
    for rows in code_rows
        code = StabilizerCode(rows)
        logical = normalize(ComplexF64[complex(i, 1 / (i+1)) for i in 1:2^code.k])
        encoded = code.encode(logical)
        @test code.decode(encoded) ≈ logical
        @test code.encode(sa.sparse(logical)) ≈ encoded
        @test all(s -> isapprox(expect(encoded, string_to_matrix(s)), 1; atol=1e-12), rows)
        for q in 1:code.k, name in ("I", "X", "Y", "Z")
            @test code.decode(code.apply(encoded, name, q)) ≈ apply(logical, Op(name, q))
        end
        measured, syndrome, correction = code.syndrome(encoded)
        @test syndrome == zeros(Int, code.m)
        @test correction.name == "I"
        @test code.decode(code.correct(measured, syndrome)) ≈ logical

        # The public binary logical construction must obey the symplectic algebra.
        _, X, Z = get_XZ_logicals(code.generator_standard, code.permute_matrix_inv)
        n = code.n
        symplectic(A, B) = mod.(A[:, 1:n] * B[:, n+1:end]' + A[:, n+1:end] * B[:, 1:n]', 2)
        @test all(iszero, symplectic(code.generator, X))
        @test all(iszero, symplectic(code.generator, Z))
        @test symplectic(X, Z) == Matrix{Int}(I, code.k, code.k)
    end

    # Every single Pauli error is correctable for the [[5,1,3]] and Steane codes.
    logical = normalize(ComplexF64[1, 0.2 + 0.7im])
    for rows in code_rows[2:3]
        code = StabilizerCode(rows; d=3)
        encoded = code.encode(logical)
        for qubit in 1:code.n, name in ("X", "Y", "Z")
            errored = apply(encoded, Op(name, qubit))
            measured, syndrome, correction = code.syndrome(errored)
            expected = [round(Int, (1 - expect(errored, string_to_matrix(s))) / 2) for s in rows]
            @test syndrome == expected
            @test correction isa Op
            @test all(isapprox.(expect(measured, "Z")[code.n+1:end], 1; atol=1e-12))
            decoded = code.decode(code.correct(measured, syndrome))
            @test isapprox(abs2(dot(logical, decoded)), 1; atol=1e-12)
        end
        @test isapprox(abs2(dot(logical, code.decode(code.correct(apply(encoded, Op("Y", 1)))))), 1; atol=1e-12)
    end

    repetition = StabilizerCode(code_rows[1])
    for q in 1:3
        @test isapprox(abs2(dot(logical, repetition.decode(repetition.correct(apply(repetition.encode(logical), Op("X", q)))))), 1; atol=1e-12)
    end
    @test_throws ArgumentError repetition.decode(apply(repetition.encode(logical), Op("X", 1)))
    @test_throws ArgumentError StabilizerCode(["X", "Z"])
    @test_throws ArgumentError StabilizerCode(["Z,I", "Z,I"])
    @test_throws ArgumentError StabilizerCode(["I,I"])
    @test_throws ArgumentError StabilizerCode(["X,I", "Z"])
    @test_throws ArgumentError StabilizerCode(["bad,I"])
    @test_throws ArgumentError StabilizerCode(String[])
    @test_throws DimensionMismatch repetition.encode(plus_state(2))
    @test_throws DimensionMismatch repetition.syndrome(zero_state(2))
    @test_throws DimensionMismatch repetition.correct(repetition.encode(), [0])
    @test_throws ArgumentError repetition.apply(repetition.encode(), "X", 2)

    supplied = Dict(("X", 1) => [("X", 1), ("X", 2), ("X", 3)])
    custom = StabilizerCode(code_rows[1], supplied)
    @test custom.logicals["X", 1] == supplied["X", 1]
    @test length(supplied) == 1
    @test length(custom.ops([Op("X", 1)])) == 3
    @test isempty(custom.ops("I", 1))
    physical = Op("X", 3)
    @test physical_ops(Op("X", 1), Dict(("X", 1) => physical)) === physical
    @test only(physical_ops([Op("X", 1)], Dict(("X", 1) => physical))) === physical

    steane = StabilizerCode(code_rows[3])
    unmatched = mod.(steane.generator[:, steane.n+1] + steane.generator[:, 2], 2)
    @test_throws ArgumentError steane.correct(steane.encode(), unmatched)
end

@testset "QEC density and MPS backends" begin
    rows = ["X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"]
    code = StabilizerCode(rows; d=3)
    logical_rho = ComplexF64[0.6 0.2+0.1im; 0.2-0.1im 0.4]
    encoded = code.encode(logical_rho)
    @test code.decode(encoded) ≈ logical_rho
    @test tr(encoded) ≈ 1
    for name in ("X", "Y", "Z")
        @test code.decode(code.apply(encoded, name, 1)) ≈ gates(name) * logical_rho * gates(name)'
    end
    noisy = 0.55 * encoded + 0.25 * apply(encoded, Op("X", 1)) + 0.2 * apply(encoded, Op("Y", 4))
    recovered = code.correct(noisy)
    @test code.decode(recovered) ≈ logical_rho
    @test code.correct(noisy) ≈ recovered
    measured, syndrome, _ = code.syndrome(noisy)
    @test code.decode(code.correct(measured, syndrome)) ≈ logical_rho
    @test code.decode(code.correct(noisy; noise=NoiseModel("depolarizing", 0.0))) ≈ logical_rho

    # Native MPS encoding and recovery agree with a small dense reference.
    logical_mps = apply(plus_state(N_MPS(1)), Op("S", 1))
    encoded_mps = code.encode(logical_mps)
    @test to_state(encoded_mps) ≈ code.encode(to_state(logical_mps))
    @test to_state(code.decode(encoded_mps)) ≈ to_state(logical_mps)
    @test to_state(code.decode(code.apply(encoded_mps, "Y", 1))) ≈ gate.Y * to_state(logical_mps)
    damaged = apply(encoded_mps, Op("Y", 2))
    measured_mps, bits, _ = code.syndrome(damaged)
    @test measured_mps isa it.MPS
    @test all(isapprox.(expect(measured_mps, "Z")[code.n+1:end], 1; atol=1e-10))
    corrected = code.correct(measured_mps, bits)
    @test isapprox(abs2(dot(to_state(logical_mps), to_state(code.decode(corrected)))), 1; atol=1e-10)
    decoded_mps = code.decode(code.correct(damaged))
    @test isapprox(abs2(dot(to_state(logical_mps), to_state(decoded_mps))), 1; atol=1e-10)
    @test_throws ArgumentError code.decode(damaged)

    # Density decoding retains mixed logical states even outside the zero-ancilla sector.
    repetition = StabilizerCode(["Z,Z,I", "I,Z,Z"])
    rho = repetition.encode(logical_rho)
    @test tr(repetition.decode(apply(rho, Op("X", 1)))) ≈ 1
    physical, logical_state = qec_state_prep(N_MPS(3), [3, 1], :zero_plus)
    @test to_state(reduce_and_relabel_qubits(physical, [3, 1])) ≈ to_state(logical_state)
end

@testset "QEC initial states" begin
    for n in (3, N_MPS(3))
        for initial in (:zero, :one, :plus, :minus, :zero_plus, [1, 0])
            physical, logical = qec_state_prep(n, [3, 1], initial)
            if !(n isa Int)
                physical, logical = to_state(physical), to_state(logical)
            end
            @test reduce_and_relabel_qubits(physical, [3, 1]) ≈ logical
        end
    end
    @test_throws ArgumentError qec_state_prep(3, [1, 1])
    @test_throws ArgumentError qec_state_prep(3, [4])
    @test_throws ArgumentError qec_state_prep(3, [1], [2])
    @test_throws DimensionMismatch qec_state_prep(3, [1], [0, 1])
    @test_throws ArgumentError reduce_and_relabel_qubits(plus_state(2), [1])

    code = StabilizerCode(["Z,Z,I", "I,Z,Z"])
    for n in (3, N_MPS(3))
        state = QECState(n, code.logical_indices, :plus; encoding_circuit=code.ops_encoding, cutoff=1e-12)
        encoded = state.is_mps ? to_state(state.encoded_state) : state.encoded_state
        @test code.decode(encoded) ≈ plus_state(1)
    end
end

@testset "QEC decoder callbacks across state representations" begin
    # Two adjacent X errors in the five-qubit repetition code have no matching
    # single-qubit recovery. A decoder for this error model can correct both.
    code = StabilizerCode(["Z,Z,I,I,I", "I,Z,Z,I,I", "I,I,Z,Z,I", "I,I,I,Z,Z"])
    recovery = [Op("X", 1; noisy=false), Op("X", 2; noisy=false)]
    bits = [0, 1, 0, 0]
    decoder = syndrome -> all(iszero, syndrome) ? Op[] : copy(recovery)
    logical = normalize(ComplexF64[1, 0.3+0.2im])

    for initial in (logical, to_MPS(logical, N_MPS(1); cutoff=0.0))
        encoded = code.encode(initial)
        damaged = apply(encoded, recovery)
        measured, syndrome, candidate = code.syndrome(damaged; decoder=decoder)
        @test syndrome == bits
        @test getproperty.(candidate, :qubit) == [1, 2]
        @test_throws ArgumentError code.correct(measured, syndrome)
        corrected = code.decode(code.correct(measured, syndrome; decoder=decoder))
        direct = code.decode(code.correct(damaged; decoder=decoder))
        @test isapprox(abs2(inner(initial, corrected)), 1; atol=1e-10)
        @test isapprox(abs2(inner(initial, direct)), 1; atol=1e-10)
        @test (corrected isa it.MPS) == (initial isa it.MPS)
    end

    logical_rho = logical * logical'
    encoded_rho = code.encode(logical_rho)
    noisy_rho = 0.4encoded_rho + 0.6apply(encoded_rho, recovery)
    @test code.decode(code.correct(noisy_rho; decoder=decoder)) ≈ logical_rho
    # This enters the explicit noisy-ancilla density branch, even at p=0.
    @test code.decode(code.correct(noisy_rho; decoder=decoder,
                                  noise=NoiseModel("depolarizing", 0.0))) ≈ logical_rho

    encoded = code.encode(logical)
    @test_throws DimensionMismatch code.correct(encoded, [0]; decoder=decoder)
    @test_throws ArgumentError code.correct(encoded, [0, 2, 0, 0]; decoder=decoder)
    for invalid in ("X", [Op("H", 1)], [Op("X", 6)], [Op("X", 1; control=2)],
                    [Op("X", gate.Z, 1)], Op[])
        @test_throws ArgumentError code.correct(encoded, bits; decoder=_ -> invalid)
    end

    # A valid zero-syndrome recovery need not be the identity. Do not discard
    # the decoder's explicit choice; syndrome consistency cannot certify it.
    flip_logical = _ -> [Op("X", q; noisy=false) for q in 1:5]
    flipped = code.decode(code.correct(code.encode(), zeros(Int, 4); decoder=flip_logical))
    @test abs2(dot(one_state(1), flipped)) ≈ 1
end
