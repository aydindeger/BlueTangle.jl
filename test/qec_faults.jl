@testset "Pauli fault propagation and logical classification" begin
    # Dense two-qubit conjugation is an independent reference, including the
    # phases discarded by the binary propagation routine.
    paulis = ("I", "X", "Y", "Z")
    pauli_matrix(ops, n) = foldl((P, op) -> Matrix(op.expand(n)) * P, ops;
                                init=Matrix{ComplexF64}(I, 2^n, 2^n))
    for name in ("I", "X", "Y", "Z", "H", "S", "SD", "CX", "CNOT", "CY", "CZ", "SWAP")
        pair = name in ("CX", "CNOT", "CY", "CZ", "SWAP")
        op = pair ? Op(name, 2, 1) : Op(name, 2)
        U = Matrix(op.expand(2))
        for p in paulis, q in paulis
            fault = [Op(p, 1), Op(q, 2)]
            propagated = propagate_pauli(fault, QuantumOps[op]; nqubits=2)
            expected = U * kron(gates(p), gates(q)) * U'
            actual = pauli_matrix(propagated, 2)
            @test abs(tr(actual' * expected)) ≈ 4
            @test all(!o.noisy && o.name != "I" for o in propagated)
        end
    end
    circuit = [Op("H", 1), Op("S", 2), Op("CX", 1, 2), Op("SD", 1), Op("SWAP", 1, 2)]
    U = pauli_matrix(circuit, 2)
    fault = QuantumOps[Op("X", 1), Op("Z", 1), Op("Y", 2)]
    actual = pauli_matrix(propagate_pauli(fault, circuit), 2)
    expected = U * pauli_matrix(fault, 2) * U'
    @test abs(tr(actual' * expected)) ≈ 4
    @test isempty(propagate_pauli([Op("X", 1), Op("X", 1)], Op[]))
    @test isempty(propagate_pauli(Op[], Op[]))
    @test only(propagate_pauli(Op("X", 3), [Op("H", 3)])).name == "Z"
    @test_throws ArgumentError propagate_pauli(Op("X", 3), circuit; nqubits=2)
    @test_throws ArgumentError propagate_pauli(Op[], Op[]; nqubits=-1)
    for invalid in (Op("T", 1), Op("MZ", 1), Op("X", 1; control=2),
                    Op("H", gate.S, 1), OpQC("depolarizing", 0.1, 1))
        @test_throws ArgumentError propagate_pauli(Op("X", 1), [invalid])
    end
    for invalid in (Op("H", 1), Op("X", 1; control=2), Op("X", gate.Z, 1), "X")
        @test_throws ArgumentError propagate_pauli([invalid], Op[]; nqubits=2)
    end

    code = StabilizerCode(["X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"]; d=3)
    encoded = code.encode(normalize(ComplexF64[1, 0.2+0.7im]))
    for q in 1:code.n, name in ("X", "Y", "Z")
        fault = Op(name, q)
        state = apply(encoded, fault)
        bits = [round(Int, (1-expect(state, string_to_matrix(s)))/2) for s in code.stabilizers]
        @test pauli_syndrome(code, fault) == bits
        @test_throws ArgumentError logical_pauli(code, fault)
        recovery = BlueTangle.identify_error(bits, code.generator)
        @test logical_pauli(code, [fault, recovery]) == (x=[0], z=[0])
    end
    for (name, x, z) in (("I", 0, 0), ("X", 1, 0), ("Y", 1, 1), ("Z", 0, 1))
        # Build logical Y as XZ, so its irrelevant global phase need not be
        # represented by the encoder's explicit phase operation.
        ops = name == "Y" ? vcat(code.ops("X", 1), code.ops("Z", 1)) : code.ops(name, 1)
        @test pauli_syndrome(code, ops) == zeros(Int, code.m)
        @test logical_pauli(code, ops) == (x=[x], z=[z])
    end
    for row in code.stabilizers
        stabilizer = [Op(String(name), q) for (q, name) in enumerate(split(row, ",")) if name != "I"]
        @test logical_pauli(code, stabilizer) == (x=[0], z=[0])
    end
    @test pauli_syndrome(code, [Op("X", 1), Op("X", 1)]) == zeros(Int, code.m)
    @test_throws ArgumentError pauli_syndrome(code, Op("X", code.n+1))
    @test_throws ArgumentError pauli_syndrome(code, Op("H", 1))
end
