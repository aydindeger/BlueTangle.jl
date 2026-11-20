# Clifford-Augmented MPS (CAMPS) walkthrough:
# 1. Exact Clifford+T compression versus plain MPS
# 2. Noisy circuit sampling and classical shadows via backend=:camps
# 3. 21-qubit Steane [[7,1,3]]^3 QEC at bond dimension 1
#
# Run from the repository root:
#   julia --project=. examples/camps_clifford_t.jl

using BlueTangle, LinearAlgebra, Random

rng = MersenneTwister(2026)

# ---------------------------------------------------------------------------
# 1. Bond-dimension compression on a nonlocal T-doped Clifford circuit (N = 10)
# ---------------------------------------------------------------------------
rng = MersenneTwister(7)
N = 10
ops = Op[]
for layer in 1:10
    for q in 1:N
        push!(ops, Op(rand(rng, ("H", "S")), q))
    end
    perm = randperm(rng, N)
    for k in 1:2:(N - 1)
        push!(ops, Op("CX", perm[k], perm[k + 1]))
    end
    push!(ops, Op("T", rand(rng, 1:N)))
end

sv = apply(zero_state(N), ops)
mps_plain = apply(zero_state(N_MPS(N)), ops; cutoff=1e-12)
st_camps = camps(ops; N=N, cutoff=1e-12, disentangle=true, nsweeps=2)

println("=== 1. Nonlocal T-doped Clifford scrambler (N = $N, 10 T gates) ===")
println("Plain MPS bond dimension:   ", maxlinkdim(mps_plain))
println("CAMPS MPS bond dimension:   ", maxlinkdim(st_camps))
println("CAMPS fidelity vs exact SV: ", fidelity(st_camps, sv))
println("Max |⟨Z⟩| error vs exact:   ", maximum(abs.(expect(st_camps, "Z") .- expect(sv, "Z"))))
@assert fidelity(st_camps, sv) ≈ 1.0
@assert maxlinkdim(st_camps) < maxlinkdim(mps_plain)

# ---------------------------------------------------------------------------
# 2. Unified backend=:camps API (noisy sampling + subsystem shadow)
# ---------------------------------------------------------------------------
println("\n=== 2. Noisy trajectories & classical shadow (backend=:camps) ===")
small_ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
nm = NoiseModel(Noise1("amplitude_damping", 0.2), Noise2("depolarizing", 0.08))
rho_exact = apply(to_rho(zero_state(3)), small_ops; noise=nm)
m_camps = measure(small_ops, 3000; backend=:camps, options=Options(noise=nm))
println("Exact noisy ⟨Z⟩: ", round.(expect(rho_exact, "Z"); digits=3))
println("CAMPS noisy ⟨Z⟩: ", round.(m_camps.expect; digits=3))

rho_13 = shadow(small_ops, 3000; qubits=[1, 3], backend=:camps, cutoff=1e-12)
rho_13_exact = partial_trace(apply(zero_state(3), small_ops), [1, 3])
println("Subsystem shadow trace error: ", round(norm(Matrix(rho_13) - Matrix(rho_13_exact)); digits=4))

# ---------------------------------------------------------------------------
# 3. 21-qubit Steane QEC (3 logical blocks) at chi = 1
# ---------------------------------------------------------------------------
println("\n=== 3. 21-qubit Steane(3) QEC with CAMPS ===")
code = steane_code(3)
decoder = steane_decoder(code)

enc = code.encode(CAMPSState(3))
enc = code.apply(enc, "H", 1)
enc = code.apply(enc, "CX", 1, 2)
enc = code.apply(enc, "CX", 2, 3)

damaged = apply(enc, [Op("X", 3), Op("Y", 11), Op("Z", 19)])
rec = qec_round(code, damaged; decoder=decoder)
logical = code.decode(rec.state)

println("Encoded physical qubits:    ", get_N(enc))
println("Bond dimension after QEC:   ", maxlinkdim(rec.state))
println("Recovered state fidelity:   ", fidelity(rec.state, enc))
println("Syndrome bits (18 checks):  ", rec.syndrome)
println("Decoded logical ⟨X₁X₂X₃⟩:   ", correlation(logical, "X,X,X", [1, 2, 3]))
@assert maxlinkdim(rec.state) == 1
@assert fidelity(rec.state, enc) ≈ 1.0
@assert correlation(logical, "X,X,X", [1, 2, 3]) ≈ 1.0
