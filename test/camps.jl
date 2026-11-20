@testset "CAMPS tableau, synthesis, and Clifford circuits" begin
    rng = MersenneTwister(2026)

    # Constructors, validation, copy, and show
    st0 = CAMPSState(4)
    @test get_N(st0) == 4
    @test maxlinkdim(st0) == 1
    @test norm(st0) ≈ 1.0
    @test entanglement_entropy(st0) ≈ 0.0
    @test occursin("CAMPSState(N=4, maxlinkdim=1)", sprint(show, st0))
    @test occursin("CliffordTableau(N=4)", sprint(show, st0.tableau))
    @test occursin("C† X_1 C", sprint(show, MIME("text/plain"), st0.tableau))
    @test_throws ArgumentError CAMPSState(0)
    @test_throws ArgumentError CliffordTableau(0)
    @test_throws DimensionMismatch CAMPSState(CliffordTableau(3), zero_state(N_MPS(4)))
    @test_throws ArgumentError camps(Op[]; N=0)
    @test_throws ArgumentError camps([Op("H", 3)]; N=2)

    # Wrapping an existing MPS preserves state without aliasing
    mps_init = apply(zero_state(N_MPS(3)), [Op("RX(0.5)", 1), Op("CX", 1, 2)])
    st_wrap = CAMPSState(mps_init)
    @test fidelity(st_wrap, mps_init) ≈ 1.0 atol=1e-12
    st_after = apply(st_wrap, Op("H", 1))
    @test fidelity(st_wrap, mps_init) ≈ 1.0 atol=1e-12 # caller not mutated
    @test !isapprox(fidelity(st_after, mps_init), 1.0; atol=1e-2)

    # Every built-in 1-qubit Clifford gate
    for gname in ("I", "X", "Y", "Z", "H", "S", "SD", "SX", "XSQRT",
                  "SXDG", "SXD", "XSQRTD", "XSQRTDG", "HSP", "HY")
        ops = [Op("RX(0.4)", 1), Op("RY(0.9)", 2), Op("CX", 1, 2), Op(gname, 1), Op(gname, 2)]
        st = camps(ops; N=2)
        sv = apply(zero_state(2), ops)
        @test same_up_to_global_phase(to_state(st), sv)[1]
        @test to_rho(st) ≈ sv * sv' atol=1e-12
        @test fidelity(st, sv) ≈ 1.0 atol=1e-12
    end

    # Every built-in 2-qubit Clifford gate in both directions
    for gname in ("CI", "CX", "CNOT", "CZ", "CY", "SWAP", "ISWAP", "FSWAP", "ECR", "H2")
        ops = [Op("RX(0.7)", 1), Op("RY(1.1)", 2), Op(gname, 1, 2), Op(gname, 2, 1)]
        st = camps(ops; N=2)
        sv = apply(zero_state(2), ops)
        @test same_up_to_global_phase(to_state(st), sv)[1]
        @test fidelity(st, sv) ≈ 1.0 atol=1e-12
    end

    # Deep random Clifford circuit stays at chi=1 and synthesizes accurately
    cliff_1q = ("H", "S", "SD", "X", "Y", "Z", "SX", "SXDG", "HSP", "HY")
    cliff_2q = ("CX", "CZ", "CY", "SWAP", "ISWAP", "FSWAP", "ECR")
    ops6 = Op[]
    for _ in 1:80
        if rand(rng, Bool)
            push!(ops6, Op(rand(rng, cliff_1q), rand(rng, 1:6)))
        else
            a = rand(rng, 1:6)
            b = rand(rng, setdiff(1:6, a))
            push!(ops6, Op(rand(rng, cliff_2q), a, b))
        end
    end
    st6 = camps(ops6; N=6)
    sv6 = apply(zero_state(6), ops6)
    @test maxlinkdim(st6) == 1
    @test entanglement_entropy(st6) ≈ 0.0 atol=1e-12
    @test same_up_to_global_phase(to_state(st6), sv6)[1]
    @test same_up_to_global_phase(to_state(to_MPS(st6; cutoff=0.0)), sv6)[1]
    @test all(op -> !op.noisy, clifford_ops(st6))

    # Parameterized gates at exact Clifford angles stay at chi=1
    exact_cliff_param = [
        Op("RX(0.5pi)", 1), Op("RY(-0.5pi)", 2), Op("RZ(pi)", 3), Op("P(0.5pi)", 4),
        Op("RZZ(0.5pi)", 1, 2), Op("RXX(0.25pi)", 2, 3), Op("RXY(-0.5pi)", 3, 4), Op("CP(pi)", 1, 4),
        Op("Xbasis", gate.H, 1), Op("Ybasis", gate.HSP, 2), Op("Zbasis", gate.I, 3),
    ]
    st_cp = camps(exact_cliff_param; N=4)
    sv_cp = apply(zero_state(4), exact_cliff_param)
    @test maxlinkdim(st_cp) == 1
    @test same_up_to_global_phase(to_state(st_cp), sv_cp)[1]

    # Large 30-qubit Clifford scrambler stays at chi=1
    ops30 = Op[]
    for _ in 1:200
        a = rand(rng, 1:30)
        b = rand(rng, setdiff(1:30, a))
        push!(ops30, Op(rand(rng, ("H", "S")), a))
        push!(ops30, Op(rand(rng, ("CX", "CZ")), a, b))
    end
    st30 = camps(ops30; N=30)
    @test maxlinkdim(st30) == 1
    @test entanglement_entropy(st30) ≈ 0.0 atol=1e-12
end

@testset "CAMPS non-Clifford gates, KAK, controlled gates, and disentangling" begin
    rng = MersenneTwister(2027)

    # 1-qubit and 2-qubit parameterized & special gates
    batches = [
        [Op("T", 1), Op("TD", 2), Op("P(0.4)", 1), Op("U1(0.7)", 2), Op("U2(0.3, 1.1)", 1), Op("U3(0.8, 0.4, -0.6)", 2)],
        [Op("RX(0.5)", 1), Op("RY(0.9)", 2), Op("RZ(1.3)", 1), Op("CS", 1, 2), Op("CSD", 2, 1)],
        [Op("RX(0.5)", 1), Op("RY(0.9)", 2), Op("CT", 1, 2), Op("CTD", 2, 1), Op("CH", 1, 2)],
        [Op("RX(0.5)", 1), Op("RY(0.9)", 2), Op("CP(0.6)", 1, 2), Op("RXX(0.35)", 1, 2), Op("RYY(0.42)", 2, 1)],
        [Op("RX(0.5)", 1), Op("RY(0.9)", 2), Op("RZZ(0.73)", 1, 2), Op("RXY(0.61)", 1, 2), Op("GIVENS(0.55)", 1, 2)],
        [Op("RX(0.5)", 1), Op("RY(0.9)", 2), Op("FSIM(0.4, 0.7)", 1, 2), Op("SYC", 2, 1), Op("SWAPA(0.33)", 1, 2)],
    ]
    for ops in batches
        st = camps(ops; N=2)
        sv = apply(zero_state(2), ops)
        @test same_up_to_global_phase(to_state(st), sv)[1]
    end

    # Arbitrary Haar-random U(2) and U(4) matrices (KAK decomposition) + degenerate 2q gates
    for _ in 1:15
        U1 = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        U4 = Matrix(qr(randn(rng, ComplexF64, 4, 4)).Q)
        ops = [Op("u1", U1, 1), Op("RY(0.7)", 2), Op("u4", U4, 1, 2), Op("u4_rev", U4, 2, 1)]
        @test same_up_to_global_phase(to_state(camps(ops; N=2)), apply(zero_state(2), ops))[1]
    end
    for U in (gate.CI, gate.CX, gate.CZ, gate.SWAP, gate.ISWAP, gate.FSWAP, gate.ECR, gate.SYC, gate.CH)
        ops = [Op("RX(0.4)", 1), Op("RY(0.9)", 2), Op("forced_kak", U, 1, 2)]
        @test same_up_to_global_phase(to_state(camps(ops; N=2)), apply(zero_state(2), ops))[1]
    end

    # Controlled 1-qubit gates and 3-qubit gates
    for ops in (
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("RZ(0.3)", 3), Op("CCZ", 1, 2, 3)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("RZ(0.3)", 3), Op("CCX", 1, 2, 3)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("RZ(0.3)", 3), Op("CCY", 1, 2, 3)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("RZ(0.3)", 3), Op("CSWAP", 1, 2, 3)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("U3(0.6, 0.2, 0.9)", 2; control=1)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("T", 2; control=1)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("X", 2; control=1)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("Y", 2; control=1)],
        [Op("RX(0.4)", 1), Op("RY(0.7)", 2), Op("Z", 2; control=1)],
    )
        st = camps(ops; N=3)
        sv = apply(zero_state(3), ops)
        @test same_up_to_global_phase(to_state(st), sv)[1]
    end

    # Variational gate (`op.type == "f"`) and multiplication syntax (`*`)
    vop = Op("vry", _RY, 1)
    st_v = apply(CAMPSState(2), vop; pars=(0.65,))
    sv_v = apply(zero_state(2), vop; pars=(0.65,))
    @test same_up_to_global_phase(to_state(st_v), sv_v)[1]
    st_mul = Op("T", 1) * ([("H", 1), ("CX", 1, 2)] * CAMPSState(2))
    sv_mul = Op("T", 1) * ([("H", 1), ("CX", 1, 2)] * zero_state(2))
    @test same_up_to_global_phase(to_state(st_mul), sv_mul)[1]

    # Free-zero Pauli reduction: 10 non-local pulled-back T gates + deep Clifford scrambler stays at chi=1!
    ops_fz = Op[Op("H", q) for q in 1:10]
    for q in 1:9
        push!(ops_fz, Op("CX", q, q + 1))
    end
    for q in 1:10
        push!(ops_fz, Op("T", q))
    end
    for layer in 1:6
        for q in 1:10
            push!(ops_fz, Op(isodd(layer + q) ? "H" : "S", q))
        end
        for q in (isodd(layer) ? 1 : 2):2:9
            push!(ops_fz, Op("CX", q, q + 1))
        end
    end
    st_fz = camps(ops_fz; N=10)
    sv_fz = apply(zero_state(10), ops_fz)
    @test maxlinkdim(st_fz) == 1
    @test fidelity(st_fz, sv_fz) ≈ 1.0 atol=1e-11

    # T-doped Clifford scrambler: CAMPS compresses dramatically vs plain MPS
    ops_td = Op[]
    for layer in 1:8
        for q in 1:8
            push!(ops_td, Op(rand(rng, ("H", "S")), q))
        end
        for q in (isodd(layer) ? 1 : 2):2:7
            push!(ops_td, Op("CX", q, q + 1))
        end
        push!(ops_td, Op("T", rand(rng, 1:8)))
        push!(ops_td, Op("T", rand(rng, 1:8)))
    end
    sv_td = apply(zero_state(8), ops_td)
    mps_td = apply(zero_state(N_MPS(8)), ops_td; cutoff=1e-12)
    st_td = camps(ops_td; N=8, cutoff=1e-12)
    st_td_nodis = camps(ops_td; N=8, cutoff=1e-12, disentangle=false)
    @test fidelity(st_td, sv_td) ≈ 1.0 atol=1e-10
    @test fidelity(st_td_nodis, sv_td) ≈ 1.0 atol=1e-10
    @test maxlinkdim(st_td) < maxlinkdim(mps_td)

    # Explicit maxdim budget & keyword validation
    st_cap = camps(ops_td; N=8, maxdim=2, cutoff=1e-12, nsweeps=2)
    @test maxlinkdim(st_cap) <= 2
    @test_throws ArgumentError apply(CAMPSState(2), Op("T", 1); maxdim=0)
    @test_throws ArgumentError apply(CAMPSState(2), Op("T", 1); cutoff=-1e-5)
    @test_throws ArgumentError apply(CAMPSState(2), Op("T", 1); nsweeps=-1)
    @test_throws Any apply(CAMPSState(2), Op("T", 1); bogus_kw=1)
end

@testset "CAMPS observables, sampling, measurements, OpF, OpQC, and noise" begin
    rng = MersenneTwister(2028)

    # Observables match exact statevector across all overloads
    ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("RY(0.6)", 3), Op("CZ", 2, 3), Op("S", 2)]
    st = camps(ops; N=3)
    sv = apply(zero_state(3), ops)

    for axis in ("I", "X", "Y", "Z", "P0", "P1")
        @test expect(st, axis) ≈ expect(sv, axis) atol=1e-12
    end
    @test expect(st, "X,Y,Z") ≈ correlation(sv, "X,Y,Z", [1, 2, 3]) atol=1e-12
    @test expect(st, "Z", 2) ≈ expect(sv, Op("Z", 2)) atol=1e-12
    @test expect(st, Op("X", 1)) ≈ expect(sv, Op("X", 1)) atol=1e-12
    @test expect(st, Op("CX", 1, 2)) ≈ expect(sv, Op("CX", 1, 2)) atol=1e-12
    @test expect(st, Op("Z", 2; control=1)) ≈ expect(sv, Op("Z", 2; control=1)) atol=1e-12
    @test expect(st, gate.Z ⊗ gate.X ⊗ gate.Y) ≈ expect(sv, gate.Z ⊗ gate.X ⊗ gate.Y) atol=1e-12
    @test correlation(st, "X,Z", [1, 3]) ≈ correlation(sv, "X,Z", [1, 3]) atol=1e-12
    @test correlation(st, "P0,P1", [1, 2]) ≈ correlation(sv, "P0,P1", [1, 2]) atol=1e-12
    @test correlation(st, [1, 2]) ≈ correlation(sv, [1, 2]) atol=1e-12

    # Exact and finite-shot sampling without mutating input state
    before = to_state(st)
    bits_exact, probs_exact = sample_exact(st)
    @test bits_exact == collect(0:7)
    @test probs_exact ≈ abs2.(sv) atol=1e-12
    @test length(sample_bit(st, 5)) == 5
    @test length(sample_bit(st, 5; basis="X")) == 5
    @test same_up_to_global_phase(to_state(st), before)[1]

    # Sample distribution matches exact Born probabilities
    m_sampled = measure(st, 2000)
    @test m_sampled.expect ≈ expect(sv, "Z") atol=0.12
    p0_est, p1_est = prob(st, 1000, 1)
    @test abs(p0_est - (1 + expect(sv, Op("Z", 1))) / 2) < 0.12
    @test p0_est + p1_est ≈ 1.0

    # Large 72-qubit stabilizer sampling (BigInt outcomes, fast stabilizer path)
    st_large = apply(CAMPSState(72), [Op("X", q) for q in 1:72])
    m_large = measure(st_large, 4)
    @test m_large.bitstr == [big(2)^72 - 1]
    @test m_large.expect == -ones(72)
    @test isempty(m_large.mid_measurements)
    @test_throws ArgumentError sample_exact(st_large)
    @test_throws ArgumentError measure(st_large)

    # Mid-circuit measurements, born_measure_Z, RESET, and ifOp (including nested branch measurements)
    st_bell = camps([Op("H", 1), Op("CX", 1, 2)]; N=2)
    collapsed, b1 = born_measure_Z(st_bell, 1)
    @test maxlinkdim(collapsed) == 1
    @test expect(collapsed, "Z") ≈ fill(1 - 2 * b1, 2) atol=1e-12

    st_reset, rec_reset = apply(CAMPSState(1), [Op("X", 1), Op("MZ", 1), Op("RES", 1), Op("MZ", 1)];
                                track_measurements=true)
    @test rec_reset == [1, 0]
    @test expect(st_reset, "Z", 1) ≈ 1.0 atol=1e-12

    cond = ifOp("MZ", 1, [Op("MX", 2), Op("MZ", 1)], [Op("X", 2), Op("MZ", 2)])
    _, rec0 = apply(CAMPSState(2), [cond]; track_measurements=true)
    @test length(rec0) == 3 && rec0[1] == 0 && rec0[3] == 0
    st1, rec1 = apply(CAMPSState(2), [Op("X", 1), cond]; track_measurements=true)
    @test rec1 == [1, 1]
    @test expect(st1, "Z") ≈ [-1.0, -1.0] atol=1e-12

    # OpF variants (block, callback, 1q/2q matrix, unsupported 3q dense matrix)
    block = Op("block", [Op("H", 1), Op("CX", 1, 2)])
    @test fidelity(apply(CAMPSState(2), [block]), apply(zero_state(2), [block])) ≈ 1.0 atol=1e-12
    cb = Op("flip1", s -> apply(s, Op("X", 1)))
    @test expect(apply(CAMPSState(2), [cb]), "Z") ≈ [-1.0, 1.0] atol=1e-12
    @test expect(apply(CAMPSState(1), Op("mat1", gate.X)), "Z") ≈ [-1.0] atol=1e-12
    @test_throws ArgumentError apply(CAMPSState(3), Op("mat3", Matrix{ComplexF64}(I, 8, 8)))

    # Noise channels: deterministic & stochastic (Pauli, coherent, amplitude & phase damping, 2q noise)
    damping1 = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    @test expect(apply(CAMPSState(1), [Op("X", 1)]; noise=damping1), "Z") ≈ [1.0] atol=1e-12
    @test expect(apply(CAMPSState(1), OpQC("bit_flip", 1.0, 1)), "Z") ≈ [-1.0] atol=1e-12
    @test expect(apply(CAMPSState(2), [Op("X", 1), Op("X", 2), OpQC("amplitude_damping", 1.0, 1, 2)]), "Z") ≈ [1.0, 1.0] atol=1e-12

    # Statistical agreement with exact density matrix under amplitude damping + depolarizing
    noisy_ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2)]
    nm_ad = NoiseModel(Noise1("amplitude_damping", 0.35), Noise2("depolarizing", 0.15))
    rho_exact = apply(to_rho(zero_state(2)), noisy_ops; noise=nm_ad)
    m_camps = measure(noisy_ops, 2000; backend=:camps, options=Options(noise=nm_ad))
    @test m_camps.expect ≈ expect(rho_exact, "Z") atol=0.08
end

@testset "CAMPS unified backend API (apply, measure, sample, shadow)" begin
    ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 3)]
    sv = apply(zero_state(3), ops)

    # apply, measure, and sample with backend=:camps and alias :clifford_mps
    for bk in (:camps, :clifford_mps)
        st_bk = apply(ops; backend=bk, cutoff=1e-12, maxdim=16, disentangle=true, nsweeps=2)
        @test st_bk isa CAMPSState
        @test fidelity(st_bk, sv) ≈ 1.0 atol=1e-11
        m = measure(ops, 1500; backend=bk, cutoff=1e-12, maxdim=16, disentangle=true, nsweeps=2)
        @test m.expect ≈ expect(sv, "Z") atol=0.10
        @test isempty(m.mid_measurements)
        draws = sample(ops, 50; backend=bk)
        @test length(draws) == 50
    end

    # layout routing, idle qubits, readout noise, twirling, measurement_basis
    layout = Layout(ones(Int, 1, 3))
    st_routed = apply(ops; backend=:camps, layout=layout)
    @test to_rho(st_routed) ≈ sv * sv' atol=1e-11

    @test measure([Op("X", 1)], 10; N=3, backend=:camps).bitstr == [4]
    @test measure(Op[], 10; N=2, backend=:camps).bitstr == [0]
    @test measure([Op("I", 1)], 20; backend=:camps,
                  options=Options(readout_noise=Noise1("bit_flip", 1.0))).bitstr == [1]
    @test measure([Op("H", 1)], 20; backend=:camps,
                  options=Options(measurement_basis="X")).bitstr == [0]
    @test measure([Op("H", 1), Op("CX", 1, 2)], 20; backend=:camps,
                  options=Options(twirl=true)).expect ≈ [0.0, 0.0] atol=0.6

    # Mid-circuit measurement records via measure(...).mid_measurements
    m_mid = measure([Op("X", 1), Op("MZ", 1), Op("RES", 1), Op("MZ", 1)], 3; backend=:camps)
    @test m_mid.mid_measurements == fill([1, 0], 3)
    @test m_mid.bitstr == [0]

    # Classical shadows on :camps (full and subsystem)
    ghz = [Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
    rho_13 = shadow(ghz, 4000; qubits=[1, 3], backend=:camps, cutoff=1e-12)
    rho_13_exact = partial_trace(apply(zero_state(3), ghz), [1, 3])
    @test size(rho_13) == (4, 4)
    @test tr(rho_13) ≈ 1.0 atol=1e-12
    @test norm(Matrix(rho_13) - Matrix(rho_13_exact)) < 0.12
end

@testset "CAMPS QEC integration (StabilizerCode and qec_round)" begin
    # 1. Single Steane block with non-Clifford logical state
    code1 = steane_code(1)
    log_ops = [Op("U3(0.7, 0.4, 0.2)", 1)]
    log_sv = apply(zero_state(1), log_ops)
    log_camps = camps(log_ops; N=1)

    enc = code1.encode(log_camps)
    @test get_N(enc) == 7
    @test maxlinkdim(enc) == 1 # Steane encoder is pure Clifford!

    # Logical Clifford gates preserve chi=1
    enc_h = code1.apply(enc, "H", 1)
    @test maxlinkdim(enc_h) == 1
    dec_h = code1.decode(enc_h)
    @test fidelity(dec_h, apply(log_sv, Op("H", 1))) ≈ 1.0 atol=1e-11

    # Single-qubit Pauli errors recovered by qec_round AND by code.syndrome / code.correct
    for err in (Op("X", 3), Op("Y", 5), Op("Z", 7))
        damaged = apply(enc, err)
        @test_throws ArgumentError code1.decode(damaged)

        round_res = qec_round(code1, damaged)
        @test fidelity(round_res.state, enc) ≈ 1.0 atol=1e-11
        @test fidelity(code1.decode(round_res.state), log_sv) ≈ 1.0 atol=1e-11

        anc_state, syn, _ = code1.syndrome(damaged)
        corr = code1.correct(anc_state, syn)
        @test fidelity(code1.decode(corr), log_sv) ≈ 1.0 atol=1e-11
    end

    # 2. Multi-block Steane(2) (14 qubits) with non-Clifford logical entanglement
    code2 = steane_code(2)
    dec2 = steane_decoder(code2)
    log2_ops = [Op("RY(0.6)", 1), Op("H", 2), Op("T", 2)]
    log2_sv = apply(zero_state(2), vcat(log2_ops, [Op("CX", 1, 2)]))
    enc2 = code2.encode(camps(log2_ops; N=2))
    enc2 = code2.apply(enc2, "CX", 1, 2)
    @test maxlinkdim(enc2) == 1
    damaged2 = apply(enc2, [Op("X", 2), Op("Z", 11)])
    rec2 = qec_round(code2, damaged2; decoder=dec2)
    @test fidelity(rec2.state, enc2) ≈ 1.0 atol=1e-10
    @test fidelity(code2.decode(rec2.state), log2_sv) ≈ 1.0 atol=1e-10

    # 3. 3-block Steane(3) (21 qubits) logical GHZ at chi=1
    code3 = steane_code(3)
    dec3 = steane_decoder(code3)
    enc3 = code3.encode(CAMPSState(3))
    enc3 = code3.apply(enc3, "H", 1)
    enc3 = code3.apply(enc3, "CX", 1, 2)
    enc3 = code3.apply(enc3, "CX", 2, 3)
    @test maxlinkdim(enc3) == 1
    damaged3 = apply(enc3, [Op("X", 3), Op("Y", 11), Op("Z", 19)])
    rec3 = qec_round(code3, damaged3; decoder=dec3)
    @test maxlinkdim(rec3.state) == 1
    @test fidelity(rec3.state, enc3) ≈ 1.0 atol=1e-12
    dec_ghz = code3.decode(rec3.state)
    @test get_N(dec_ghz) == 3
    @test maxlinkdim(dec_ghz) == 1
    @test correlation(dec_ghz, "X,X,X", [1, 2, 3]) ≈ 1.0 atol=1e-12
    @test correlation(dec_ghz, "Z,Z", [1, 2]) ≈ 1.0 atol=1e-12
    @test correlation(dec_ghz, "Z,Z", [2, 3]) ≈ 1.0 atol=1e-12
end
