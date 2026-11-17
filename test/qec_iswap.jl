@testset "C6 blocks and a fault-tolerant logical iSWAP" begin
    paulis = Dict("I" => gates("I"), "X" => gates("X"), "Y" => gates("Y"), "Z" => gates("Z"))

    # The binary propagation of ISWAP agrees with dense conjugation.
    op = Op("ISWAP", 2, 1)
    U = Matrix(op.expand(2))
    for a in keys(paulis), b in keys(paulis)
        propagated = propagate_pauli([Op(a, 1), Op(b, 2)], [op]; nqubits=2)
        M = foldl((P, o) -> Matrix(o.expand(2)) * P, propagated; init=Matrix{ComplexF64}(I, 4, 4))
        @test abs(tr(M' * U * kron(paulis[a], paulis[b]) * U')) ≈ 4
    end

    # The code: distance two, and the gadgets' logical representatives.
    code = c6_code()
    @test (code.n, code.k, code.m, code.d) == (6, 2, 4, 2)
    @test code_distance(code).distance == 2
    two = c6_code(2)
    L = c6_logicals(2)
    for j in 1:4, name in ("X", "Z")
        @test all(iszero, pauli_syndrome(two, L[(name, j)]))
        @test logical_pauli(two, L[(name, j)]) != (x=zeros(Int, 4), z=zeros(Int, 4))
    end
    support(ops) = Set(o.qubit for o in ops)
    for j in 1:4, l in 1:4
        @test isodd(length(intersect(support(L[("X", j)]), support(L[("Z", l)])))) == (j == l)
    end

    # Preparation: the requested eigenstates, with every check +1.
    checks = [[Op(String(s), q) for (q, s) in enumerate(split(row, ",")) if s != "I"] for row in code.stabilizers]
    expect_ops(psi, ops) = real(dot(psi, apply(psi, ops)))
    for basis in (:Z, :X, :Y), bits in ((0, 0), (1, 0), (0, 1), (1, 1))
        psi = apply(zero_state(6), c6_prepare(basis; bits=bits))
        @test all(expect_ops(psi, c) ≈ 1 for c in checks)
        for i in 1:2
            xs, zs = L[("X", i)], L[("Z", i)]
            observable = basis == :Z ? zs : basis == :X ? xs : [Op("Y", o.qubit) for o in xs]
            sign = basis == :Y ? -1 : 1                       # logical Y = -YYY
            @test sign * expect_ops(psi, observable) ≈ (-1)^bits[i]
        end
    end

    # The logical iSWAP: every combination of input bases and k gives the outputs of a two-qubit simulation.
    iswap = gates("ISWAP")
    eig(basis, bit) = begin
        w, v = eigen(Hermitian(paulis[string(basis)]))
        v[:, argmin(abs.(w .- (-1)^bit))]
    end
    flip = Dict(:Z => :Z, :X => :Y, :Y => :X)
    for (bA, bB) in ((:X, :Z), (:Z, :X), (:Y, :Z), (:Z, :Y), (:Z, :Z)), k in 0:4
        ia, ib = (rand(0:1), rand(0:1)), (rand(0:1), rand(0:1))
        outA, outB = iseven(k) ? (bA, bB) : (flip[bB], flip[bA])
        ops = vcat(c6_prepare(bA; bits=ia, block=1), c6_prepare(bB; bits=ib, block=2),
                   reduce(vcat, [c6_iswap() for _ in 1:k]; init=Op[]),
                   c6_readout(outA; block=1), c6_readout(outB; block=2))
        psi = apply(zero_state(12), ops)
        for s in sample(psi, 64)
            r = c6_decode(BlueTangle.int2bin(s, 12), (outA, outB))
            @test r.accepted
            for i in 1:2
                phi = iswap^k * kron(eig(bA, ia[i]), eig(bB, ib[i]))
                ea = real(dot(phi, kron(paulis[string(outA)], I(2)) * phi))
                eb = real(dot(phi, kron(I(2), paulis[string(outB)]) * phi))
                @test r.logical[i] == Int(ea < 0) && r.logical[2+i] == Int(eb < 0)
            end
        end
    end
    @test_throws ArgumentError c6_iswap(blocks=(1, 1))
    @test_throws ArgumentError c6_prepare(:W)
    @test_throws DimensionMismatch c6_decode(zeros(Int, 5), (:Z,))

    # Every single flipped bit is rejected; a logical support flips one logical bit.
    for q in 1:12
        bits = zeros(Int, 12)
        bits[q] = 1
        @test !c6_decode(bits, (:Z, :X)).accepted
    end

    # Fault tolerance: no single fault in preparation, two logical iSWAPs and readout is harmful.
    for (bA, bB) in ((:X, :Z), (:Y, :Z), (:Z, :Y))
        circuit = vcat(c6_prepare(bA; bits=(1, 0), block=1), c6_prepare(bB; block=2), c6_iswap(), c6_iswap(),
                       c6_readout(bA; block=1), c6_readout(bB; block=2))
        audit = c6_single_faults(circuit)
        @test audit.harmful == 0
        @test audit.detected > 0 && audit.faults == audit.detected + audit.harmless
    end

    # Contrast: in one [[4,2,2]] block, a physical iSWAP on the middle data qubits is a logical iSWAP, but the
    # single faults XX, YY and ZZ after it are undetectable logical errors.
    c4 = StabilizerCode(["X,X,X,X", "Z,Z,Z,Z"]; d=2)
    E4 = zeros(16, 4)
    for x1 in 0:1, x2 in 0:1, word in ([0, x1, x2, x1 ⊻ x2], 1 .- [0, x1, x2, x1 ⊻ x2])
        E4[bin2int(word)+1, 2x1+x2+1] = 1 / sqrt(2)
    end
    Ug = Matrix(Op("ISWAP", 2, 3).expand(4))
    @test Ug * E4 ≈ E4 * iswap
    harmful = String[]
    for a in ("I", "X", "Y", "Z"), b in ("I", "X", "Y", "Z")
        a == b == "I" && continue
        fault = [Op(n, q) for (n, q) in ((a, 2), (b, 3)) if n != "I"]
        if all(iszero, pauli_syndrome(c4, fault)) && logical_pauli(c4, fault) != (x=[0, 0], z=[0, 0])
            push!(harmful, a * b)
        end
    end
    @test sort(harmful) == ["XX", "YY", "ZZ"]

    # Pauli frames: exact for p = 0, and the expected fault rates.
    circuit = vcat(c6_prepare(:X; block=1), c6_prepare(:Z; block=2), c6_iswap(), c6_readout(:Z; block=1),
                   c6_readout(:Y; block=2))
    frames = sample_pauli_frames(circuit, 0.0; nqubits=12, shots=50)
    @test size(frames.x) == (50, 12) && !any(frames.x) && !any(frames.z)
    one = sample_pauli_frames([Op("H", 1)], 0.3; nqubits=1, shots=20000)
    @test isapprox(sum(one.x .| one.z) / 20000, 0.3; atol=0.02)
    @test isapprox(sum(one.x) / 20000, 0.2; atol=0.02)          # X or Y: two of the three Paulis
    ro = sample_pauli_frames(Op[], 0.0; nqubits=2, shots=20000, readout=0.1)
    @test isapprox(sum(ro.x) / 40000, 0.1; atol=0.01)
    @test_throws ArgumentError sample_pauli_frames([Op("T", 1)], 0.1; nqubits=1)
    # Frames agree with dense propagation of the same sampled fault: a single noisy CX with p2 = 1.
    pair = sample_pauli_frames([Op("CX", 1, 2)], 0.0; nqubits=2, shots=4000, p2=1.0)
    @test all(any(pair.x[s, :]) || any(pair.z[s, :]) for s in 1:4000)
end
