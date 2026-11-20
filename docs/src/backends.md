# Backends

BlueTangle separates **circuit description** (`Vector{QuantumOps}`, tuple
shorthands such as `[("H", 1), ("CX", 1, 2)]`, or a mix of both) from
**simulation representation**. The same operation list can be executed on five
backends—either by passing an initial state (`zero_state(N; backend=...)`) to
`apply(state, ops)`, or by passing `backend=...` directly to `apply(ops; backend=...)`,
`measure`, `sample`, and `shadow`:

- `:statevector` (aliases `:state_vector`, `:state`, `:sv`)
- `:density_matrix` (aliases `:density`, `:dm`, `:rho`)
- `:tensor` (aliases `:mps`, `:tensor_network`, `:tn`)
- `:camps` (aliases `:clifford_mps`, `:clifford-mps`)
- `:phasepoly` (aliases `:phase_poly`, `:phase-poly`; optional extension)

## Overview and comparison

| Backend | State representation | Typical scale | Noise & mid-circuit measurements | Gate set | Backend keywords |
|:--|:--|:--|:--|:--|:--|
| `:statevector` | `SparseVector{ComplexF64}` (`2^N`) | `N ≲ 18–20` | Single stochastic trajectory per execution | Universal (all `Op`, `ifOp`, `OpQC`, `OpF`) | `pars` (variational) |
| `:density_matrix` | `SparseMatrixCSC{ComplexF64}` (`2^N × 2^N`) | `N ≲ 10–12` | Exact Kraus sum for noise; `apply` averages measurements, `track_measurements=true` (and `measure(...).mid_measurements`) conditions on sampled bits | Universal (all `Op`, `ifOp`, `OpQC`, `OpF`) | `pars` (variational) |
| `:tensor` (`:mps`) | `ITensorMPS.MPS` (`O(N χ²)`) | Large `N` if `χ` stays small (1D area-law) | Single stochastic trajectory per execution | Universal (all `Op`, `ifOp`, `OpQC`, `OpF`) | `cutoff`, `maxdim`, `pars` |
| `:camps` | `CAMPSState`: `CliffordTableau` (`O(N²)`) + `ITensorMPS.MPS` (`O(N χ²)`) | Large `N` for Clifford-rich, QEC, and moderate-`T` circuits | Single stochastic trajectory per execution (Pauli noise at zero MPS cost) | Universal 1q/2q (`U(2)`, `U(4)` KAK), controlled 1q, `CCX`/`CCZ`/`CCY`/`CSWAP`, `ifOp`, `OpQC`, `OpF` (≤2q dense) | `cutoff`, `maxdim`, `disentangle`, `nsweeps`, `pars` |
| `:phasepoly` | `PhasePolyState`: affine permutations + phase polynomials + Gauss sums | Large `N` (`N ≥ 80`) in single-branch noiseless regime | **Noiseless & measurement-free only** (raises `ArgumentError` on noise/measurement) | Restricted: `H`, `Y`, `HSP`, `HY`, `X`, `CX`, `SWAP`, diagonal Clifford/`T`/`√T`, `P`/`U1`/`RZ` (`kπ/8`), `CP` (`kπ/4`), `CCZ`, `CCS`, `CCSdg` | `max_h`, `max_cccz`, `max_ccs`, `max_nonzeros` |

All built-in backends share BlueTangle's conventions: qubit indices start at `1`,
qubit `1` is the most significant bit (MSB) in statevectors and integer outcomes,
and operations execute in their exact listed order unless `route(ops, layout)` or
`layout=...` is requested. Operation lists can be written either with explicit
`Op(...)` constructors, with tuple shorthands such as `[("H", 1), ("CX", 1, 2)]`
(where each tuple `t` expands via `Op(t...)`), or as a mix of both (e.g.
combining tuples with `ifOp(...)` or `Op("X", 2; control=1)`).

The example below executes the same 3-qubit circuit on all four built-in
backends—using tuple shorthand directly—and verifies that their observables agree:

```@example backends_overview
using BlueTangle

ops = [("H", 1), ("T", 1), ("CX", 1, 2), ("RY(0.7)", 3), ("CZ", 2, 3)]

sv   = apply(zero_state(3), ops)
rho  = apply(to_rho(zero_state(3)), ops)
mps  = apply(zero_state(N_MPS(3)), ops; cutoff=1e-12)
st_c = apply(CAMPSState(3), ops; cutoff=1e-12)

for s in (rho, mps, st_c)
    @assert isapprox(expect(s, "Z"), expect(sv, "Z"); atol=1e-12)
    @assert isapprox(correlation(s, "X,Z", [1, 3]), correlation(sv, "X,Z", [1, 3]); atol=1e-12)
end

results = Dict(
    bk => round.(measure(ops, 1500; backend=bk).expect; digits=3)
    for bk in (:statevector, :density_matrix, :tensor, :camps)
)
results
```

## [Statevector (`:statevector`)](@id backend-statevector)

The default backend evolves a normalized sparse statevector of length `2^N`
(`zero_state(N)`). Every gate acts via its exact sparse operator on the full
Hilbert space, preserving both relative and **global phase** to machine
precision.

When noise (`NoiseModel`, `OpQC`) or mid-circuit measurements (`MX`, `MY`, `MZ`,
`MR`, `ifOp`) are present, each `apply` call follows **one stochastic
trajectory** sampled from the Born or Kraus probabilities. Calling
`measure(ops, shots; backend=:statevector, options=Options(noise=nm))` repeats
fresh independent trajectories for every shot and stores per-shot mid-circuit
outcomes in `.mid_measurements`.

```@example backends_sv
using BlueTangle

ops = [Op("H", 1), Op("CX", 1, 2), Op("T", 2)]
sv = apply(ops; backend=:statevector)

# Exact amplitudes, expectations, and exact probabilities (omit shots)
exact_m = measure(sv)
@assert exact_m.bitstr == [0, 3]
@assert exact_m.sample ≈ [0.5, 0.5]

# One stochastic trajectory with explicit mid-circuit record
final_sv, bits = apply(zero_state(2), [ops; Op("MZ", 1); Op("MZ", 2)]; track_measurements=true)
@assert bits[1] == bits[2]
(amplitudes = Vector(sv), record = bits)
```

### Strengths and limitations (`:statevector`)

- **Strengths:** Exact amplitudes and global phase, zero truncation error,
  supports every operation (`Op`, arbitrary multi-qubit `OpF` matrices,
  variational gradients, QEC, Hamiltonian evolution), and fast for `N ≲ 16`.
- **Limitations:**
  - Memory and runtime scale as `O(2^N)` (~16 MB at `N = 20` for a dense
    `ComplexF64` vector; intractable beyond ~24–28 qubits).
  - Noisy circuits require many independent shots to converge ensemble averages
    (Monte Carlo error scales as `1 / √shots`); for `N ≲ 10`, `:density_matrix`
    gives exact noise averages with zero sampling noise.

## [Density matrix (`:density_matrix`)](@id backend-density-matrix)

The density-matrix backend evolves a `2^N × 2^N` sparse density matrix
(`to_rho(zero_state(N))`). Unlike trajectory backends, it evolves the complete
mixed-state ensemble deterministically:

- `apply(rho, ops; noise=nm)` applies exact Kraus sums `∑ K_m ρ K_m†` for both
  noise channels and **unobserved** mid-circuit measurements.
- `apply(rho, ops; noise=nm, track_measurements=true)` samples and **conditions**
  on explicit measurement outcomes (returning `(conditional_rho, bits)`), while
  still summing noise and reset Kraus branches deterministically.
- `measure(ops, shots; backend=:density_matrix, options=Options(noise=nm))`
  computes the final mixed state **once** and draws all `shots` from it when no
  explicit mid-circuit measurements, twirling, random `"R"` basis, or callbacks
  are present (otherwise it samples per-shot records into `.mid_measurements`).

```@example backends_dm
using BlueTangle, LinearAlgebra

nm = NoiseModel(Noise1("amplitude_damping", 0.3), Noise2("depolarizing", 0.1))
ops = [Op("H", 1), Op("CX", 1, 2)]

rho = apply(ops; backend=:density_matrix, noise=nm)
exact_m = measure(rho) # exact mixed-state probabilities, zero shot noise

# Compare deterministic apply vs conditioned apply(...; track_measurements=true)
rho_unobs = apply(to_rho(zero_state(2)), [ops; Op("MZ", 1)])
cond_rho, bits = apply(to_rho(zero_state(2)), [ops; Op("MZ", 1)]; track_measurements=true)

@assert tr(rho) ≈ 1
@assert real(tr(rho_unobs^2)) ≈ 0.5   # maximally mixed over {00, 11}
@assert real(tr(cond_rho^2)) ≈ 1.0    # pure conditioned on observed bit
(
    noisy_Z = round.(expect(rho, "Z"); digits=4),
    noisy_purity = round(real(tr(rho^2)); digits=4),
    conditioned_bit = bits,
)
```

### Strengths and limitations (`:density_matrix`)

- **Strengths:** Exact, shot-noise-free noise channels, partial traces, purities,
  exact readout-mitigation benchmarks, and deterministic ensemble `qec_round`
  (returning `syndrome=nothing`).
- **Limitations:**
  - Memory scales as `O(4^N)` and gate application as `O(4^N)`, limiting
    practical simulations to roughly `N ≲ 10–12` qubits.
  - Not supported by `shadow` (classical shadows rely on individual random-basis
    shot trajectories).
  - Unobserved `apply(rho, Op("MZ", q))` averages outcomes rather than keeping a
    single branch; pass `track_measurements=true` when conditioning on a
    measurement record.

## [Matrix Product States (`:tensor` / `:mps`)](@id backend-mps)

The MPS backend represents a pure state as a chain of rank-3 tensors
(`ITensorMPS.MPS`) on qubit sites created with `N_MPS(N)`. Memory scales as
`O(N χ²)` and two-qubit gates scale as `O(χ³)`, where `χ = maxlinkdim(psi)` is
the maximum bond dimension across cuts of the 1D chain.

Truncation is controlled by `cutoff` (discarded squared singular-value weight)
and `maxdim` (hard upper bound on `χ`), passed to `apply`, `measure`, `sample`,
and `shadow` (and propagated into conditional `ifOp` branches).

```@example backends_mps
using BlueTangle

# 30-qubit 1D brickwall circuit (2^30 amplitudes — infeasible for :statevector)
N = 30
ops = Op[]
for q in 1:N
    push!(ops, Op("RY(0.35)", q))
end
for q in 1:2:(N - 1)
    push!(ops, Op("CX", q, q + 1))
end
for q in 2:2:(N - 1)
    push!(ops, Op("CZ", q, q + 1))
end

psi = apply(zero_state(N_MPS(N)), ops; cutoff=1e-12, maxdim=32)
@assert get_N(psi) == 30
@assert maxlinkdim(psi) <= 4

m = measure(psi, 500)
(
    qubits = get_N(psi),
    bond_dim = maxlinkdim(psi),
    mid_cut_entropy = round(entanglement_entropy(psi); digits=4),
    Z15 = round(expect(psi, "Z", 15); digits=4),
    sampled_Z15 = round(m.expect[15]; digits=3),
)
```

Classical shadows can reconstruct reduced density matrices directly in MPS form
without ever forming a `2^N` statevector:

```@example backends_mps
ghz_ops = [Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3), Op("CX", 3, 4)]
rho_14 = shadow(ghz_ops, 2000; qubits=[1, 4], backend=:mps, cutoff=1e-12, maxdim=16)
@assert size(rho_14) == (4, 4)
round.(real.(Matrix(rho_14)); digits=2)
```

### Strengths and limitations (`:tensor` / `:mps`)

- **Strengths:** Scales to tens or hundreds of qubits for 1D low-entanglement
  (area-law) circuits, shallow circuits, variational ansätze, and subsystem
  classical shadows. Supports full noise trajectories, mid-circuit measurements,
  resets, `ifOp`, and exact `qec_round`.
- **Limitations:**
  - **Entanglement barrier:** Even pure Clifford circuits, nonlocal 2D/all-to-all
    gates, and multi-block stabilizer QEC codes can drive `χ` exponentially toward
    `2^{N/2}`. Use `:camps` when entanglement is dominated by Clifford gates.
  - **Truncation bias:** Tight `maxdim` or loose `cutoff` silently alters the
    simulated state and outcome probabilities—separate truncation error from
    finite-shot sampling noise by checking convergence in `maxdim`.
  - **Exact enumeration cap:** `sample_exact(psi)` and `measure(psi)` (without
    `shots`) materialize `2^N` amplitudes and default to `max_qubits = 16`. Pass
    a shot count (`measure(psi, shots)`) for larger registers.
  - **Single trajectory:** Like `:statevector`, each `apply` call samples one
    noise/measurement trajectory.

## [Clifford-Augmented MPS (`:camps`)](@id backend-camps)

Matrix product states struggle with volume-law stabilizer entanglement—even a
pure Clifford circuit or multi-block stabilizer code causes `:tensor` bond
dimension to grow exponentially. **Clifford-Augmented Matrix Product States
(`CAMPSState`, `backend=:camps`)** separate stabilizer entanglement from
non-Clifford ("magic") entanglement by writing an `N`-qubit pure state as

```math
|\Psi\rangle = C |\phi\rangle,
```

where:

- `state.tableau::CliffordTableau` tracks the `N`-qubit Clifford frame `C`
  (up to global phase) in `O(N²)` bits via `px[j] = C† Xⱼ C` and
  `pz[j] = C† Zⱼ C`, and
- `state.mps::it.MPS` stores the disentangled non-Clifford MPS core `|\phi\rangle`.

### How CAMPS works

1. **Clifford gates (`O(N)` cost, zero MPS work):** Every 1-qubit and 2-qubit
   Clifford gate updates `C ← U C` directly on `state.tableau` without touching
   `state.mps` or increasing its bond dimension. Parameterized gates at exact
   multiples of `π/2` (such as `RX(pi/2)`, `RZZ(pi/2)`, `RXX(pi/4)`, `CP(pi)`)
   and Pauli noise channels take this same zero-MPS-cost path.
2. **Pauli pullback for non-Clifford rotations:** A rotation `exp(-i (θ/2) P)`
   acting on `C |\phi\rangle` equals `C exp(-i (θ/2) P̃) |\phi\rangle` with
   `P̃ = C† P C` computed in `O(N)` bitvector operations.
3. **Free-zero Pauli reduction (exact, zero bond growth):** Before applying
   `exp(-i (θ/2) P̃)` to `|\phi\rangle`, CAMPS finds untouched product `|0⟩`
   sites in `|\phi\rangle`. If `P̃` has `X` or `Y` on a free zero `q₀`,
   controlled Cliffords with control `q₀` leave `|\phi\rangle` identically
   unchanged while collapsing `P̃` to a weight-1 operator `±X_{q₀}` so the
   rotation acts locally with zero bond growth. Projective Pauli measurements
   on free zeros similarly sample in `O(N²)` tableau operations and reset `q₀`
   to `|0⟩`, recovering exact Gottesman–Knill evolution at `maxlinkdim == 1`.
4. **Two-qubit Clifford disentangling sweeps:** When a multi-qubit `P̃` must be
   applied via an MPS direct sum, CAMPS sweeps across adjacent bonds `(b, b+1)`
   testing the 19 left cosets of `C₁^{⊗2} \ C₂` (identity, 9 CNOT-class, and 9
   double-CNOT/iSWAP-class cosets). Whenever a two-qubit Clifford `D` lowers
   bond von Neumann entropy, CAMPS updates `|\phi\rangle ← D |\phi\rangle` and
   `C ← C D†`, leaving `C |\phi\rangle` invariant while compressing `|\phi\rangle`.

### Supported operations and truncation controls

`CAMPSState` supports all built-in 1- and 2-qubit Cliffords (`I, X, Y, Z, H, S,
SD, SX/XSQRT, SXDG/SXD/XSQRTD/XSQRTDG, HSP, HY, H2, CI, CX/CNOT, CZ, CY, SWAP,
ISWAP, FSWAP, ECR`), all parameterized 1- and 2-qubit gates (`T, TD, P, U1, U2,
U3, RX, RY, RZ, CS, CSD, CT, CTD, CH, CP, RXX, RYY, RZZ, RXY, GIVENS, FSIM,
SYC, SWAPA`), arbitrary custom `U(2)` (via ZYZ) and `U(4)` (via exact
Kraus–Cirac magic-basis KAK) unitaries, arbitrary controlled 1q unitaries,
`CCZ`/`CCX`/`CCY`/`CSWAP`, variational `Op(name, f, q)`, mid-circuit
measurements (`MX, MY, MZ, MR`, `born_measure_Z`), `RES`/`RESET`, `ifOp`,
`OpQC`, `NoiseModel` (Pauli fast path + exact Pauli-basis Kraus sampling for
`amplitude_damping` / `phase_damping`), and `OpF` blocks/callbacks.

| Keyword | Default | Role |
|:--|:--|:--|
| `cutoff` | `1e-12` | Singular-value truncation threshold for MPS direct sums and sweeps. |
| `maxdim` | `nothing` | Optional hard cap on the core MPS bond dimension (`nothing` = uncapped). When set, direct sums temporarily retain up to `2 * maxdim` before disentangling back to `maxdim`. |
| `disentangle` | `true` | Enable free-zero Pauli reduction and two-qubit Clifford disentangling sweeps. |
| `nsweeps` | `2` | Maximum left/right two-qubit disentangling sweeps per multi-qubit non-Clifford gate or projection. |

### Example 1: Compressing a nonlocal Clifford+T circuit vs. `:tensor`

On a 10-qubit nonlocal Clifford circuit doped with 10 `T` gates, plain `:tensor`
MPS reaches the volume-law ceiling `χ = 32`, whereas `camps` stays at `χ = 1`
with unit fidelity:

```@example backends_camps
using BlueTangle, Random

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

sv        = apply(zero_state(N), ops)
mps_plain = apply(zero_state(N_MPS(N)), ops; cutoff=1e-12)
st_camps  = camps(ops; N=N, cutoff=1e-12, disentangle=true, nsweeps=2)
st_nodis  = camps(ops; N=N, cutoff=1e-12, disentangle=false)

@assert fidelity(st_camps, sv) ≈ 1
@assert isapprox(expect(st_camps, "Z"), expect(sv, "Z"); atol=1e-12)
@assert maxlinkdim(st_camps) == 1
(
    plain_mps_chi = maxlinkdim(mps_plain),
    camps_nodis_chi = maxlinkdim(st_nodis),
    camps_chi = maxlinkdim(st_camps),
    fidelity = fidelity(st_camps, sv),
)
```

When extra non-Clifford gates push past `N`, magic entanglement grows in
`state.mps`; compare two `maxdim` caps to verify convergence:

```@example backends_camps
more_ops = copy(ops)
for layer in 1:3
    for q in 1:N
        push!(more_ops, Op(rand(rng, ("H", "S")), q))
    end
    perm = randperm(rng, N)
    for k in 1:2:(N - 1)
        push!(more_ops, Op("CX", perm[k], perm[k + 1]))
    end
    push!(more_ops, Op("T", rand(rng, 1:N)))
    push!(more_ops, Op("T", rand(rng, 1:N)))
end
sv_more = apply(zero_state(N), more_ops)
st_chi2 = camps(more_ops; N=N, maxdim=2, cutoff=1e-12)
st_chi8 = camps(more_ops; N=N, maxdim=8, cutoff=1e-12)
@assert maxlinkdim(st_chi2) <= 2 && maxlinkdim(st_chi8) <= 8
(chi2_fid = round(fidelity(st_chi2, sv_more); digits=5),
 chi8_fid = round(fidelity(st_chi8, sv_more); digits=5))
```

### Example 2: Noisy trajectories, feedback, and classical shadows

```@example backends_camps_api
using BlueTangle, LinearAlgebra

ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
nm = NoiseModel(Noise1("amplitude_damping", 0.25), Noise2("depolarizing", 0.1))

rho_exact = apply(to_rho(zero_state(3)), ops; noise=nm)
m_camps   = measure(ops, 3000; backend=:camps, options=Options(noise=nm), cutoff=1e-12, maxdim=16)
rho_13    = shadow(ops, 3000; qubits=[1, 3], backend=:camps, cutoff=1e-12)

# Dynamic circuit with conditional branch, reset, and mid-circuit records
dyn_ops = [Op("H", 1), Op("MZ", 1; if0=Op("I", 2), if1=Op("X", 2)), Op("RES", 1), Op("MZ", 1)]
_, bits = apply(dyn_ops; backend=:camps, track_measurements=true)
m_dyn   = measure(dyn_ops, 5; backend=:camps)
@assert bits[2] == 0
@assert all(r -> r[2] == 0, m_dyn.mid_measurements)
@assert isapprox(m_camps.expect, expect(rho_exact, "Z"); atol=0.08)
@assert size(rho_13) == (4, 4)
(camps_noisy_Z = round.(m_camps.expect; digits=3), exact_noisy_Z = round.(expect(rho_exact, "Z"); digits=3))
```

### Example 3: 21-qubit Steane QEC at bond dimension 1

Because stabilizer encoders, transversal Clifford gates, Pauli errors, and ideal
stabilizer projections are Clifford operations, multi-block QEC simulations stay
at `maxlinkdim == 1` on `CAMPSState`:

```@example backends_camps_qec
using BlueTangle

code = steane_code(3) # 21 physical qubits, 3 logical qubits
decoder = steane_decoder(code)

enc = code.encode(CAMPSState(3))
enc = code.apply(enc, "H", 1)
enc = code.apply(enc, "CX", 1, 2)
enc = code.apply(enc, "CX", 2, 3)

damaged   = apply(enc, [Op("X", 3), Op("Y", 11), Op("Z", 19)])
recovered = qec_round(code, damaged; decoder=decoder)
logical   = code.decode(recovered.state)

@assert get_N(enc) == 21
@assert maxlinkdim(recovered.state) == 1
@assert fidelity(recovered.state, enc) ≈ 1
@assert correlation(logical, "X,X,X", [1, 2, 3]) ≈ 1
(physical_qubits = get_N(enc), bond_dim = maxlinkdim(recovered.state),
 syndrome = recovered.syndrome, logical_XXX = correlation(logical, "X,X,X", [1, 2, 3]))
```

### Inspecting `clifford_ops` and converting states

Use [`clifford_ops`](@ref) to synthesize an `O(N²)` Clifford circuit for
`state.tableau`, and `to_state`, `to_rho`, or `to_MPS` to materialize the state:

```@example backends_camps_qec
ops2 = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("S", 2)]
st2  = camps(ops2; N=2)
sv2  = apply(zero_state(2), ops2)

c_ops = clifford_ops(st2)
@assert same_up_to_global_phase(to_state(st2), sv2)[1]
@assert to_rho(st2) ≈ sv2 * sv2'
@assert abs2(inner(to_MPS(st2; cutoff=1e-12), sv2)) ≈ 1
(synthesized_clifford_gates = length(c_ops), mps_core_chi = maxlinkdim(st2))
```

!!! note "Global phase semantics"
    `CliffordTableau` tracks `C` up to a single global phase `e^{iα}`.
    Consequently, `to_state(state)` and `to_MPS(state)` equal the exact state
    up to a single global phase (`same_up_to_global_phase`), whereas `to_rho`,
    `expect`, `correlation`, `fidelity`, `sample`, `sample_bit`, `sample_exact`,
    `measure`, `prob`, `shadow`, and `qec_round` are exact and global-phase
    independent.

### Strengths and limitations (`:camps`)

- **Strengths:** Absorbs Clifford/stabilizer entanglement at zero MPS cost
  (turning volume-law `χ = 2^{N/2}` Clifford and QEC states into `χ = 1`),
  supports universal 1- and 2-qubit gates, dynamic measurements/resets, general
  Kraus noise, multi-block QEC, and a zero-ITensor fast stabilizer sampler when
  `|\phi\rangle = |0…0⟩`. See also the runnable script
  [`examples/camps_clifford_t.jl`](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/camps_clifford_t.jl).
- **Limitations:**
  - **Magic barrier:** CAMPS removes stabilizer entanglement, not non-stabilizer
    ("magic") entanglement. Once the number of independent non-Clifford gates
    substantially exceeds `N`, `state.mps` bond dimension still grows
    exponentially unless capped with `maxdim` (which introduces truncation error).
  - **Exact enumeration cap:** `sample_exact(st)` and `measure(st)` materialize
    `2^N` amplitudes and default to `max_qubits = 16`; pass a shot count for
    larger registers.
  - **Re-absorbing `C`:** `to_MPS(st)` applies `clifford_ops(st)` back onto
    `state.mps`, which can dramatically increase bond dimension if `C` is deeply
    entangling. Evaluate `expect`, `correlation`, `sample`, and `qec_round`
    directly on `CAMPSState`.
  - **Dense `OpF` matrices:** Full-register matrix `OpF(name, mat)` above 2
    qubits is unsupported on `CAMPSState`; decompose into 1- and 2-qubit gates.

### [CAMPS API reference](@id CAMPS-API-reference)

```@docs
CAMPSState
CliffordTableau
camps
clifford_ops
```

`CAMPSState` also integrates with `apply`, `measure`, `sample`, `shadow`,
`expect`, `correlation`, `fidelity`, `born_measure_Z`, `to_state`, `to_rho`,
`to_MPS` ([Functions](functions.md)), and `StabilizerCode` / `qec_round`
([Quantum Error Correction](qec.md)).

### CAMPS references

- X. Qian, J. Huang, and M. Qin, *"Augmenting Density Matrix Renormalization Group with Clifford Circuits"*, Phys. Rev. Lett. **133**, 190402 (2024), [arXiv:2405.09217](https://arxiv.org/abs/2405.09217).
- S. Masot-Llima and A. García-Sáez, *"Stabilizer Tensor Networks: Universal Quantum Simulator on a Basis of Stabilizer States"*, Phys. Rev. Lett. **133**, 230601 (2024), [arXiv:2403.08724](https://arxiv.org/abs/2403.08724).
- G. Lami, T. Haug, and J. De Nardis, *"Quantum State Designs with Clifford Enhanced Matrix Product States"*, PRX Quantum **6**, 010345 (2025), [arXiv:2404.18751](https://arxiv.org/abs/2404.18751).
- G. E. Fux, B. Béri, R. Fazio, and E. Tirrito, *"Disentangling unitary dynamics with classically simulable quantum circuits"*, [arXiv:2410.09001](https://arxiv.org/abs/2410.09001) (2024).
- Z. Liu and B. K. Clark, *"Classical simulability of Clifford+T circuits with Clifford-augmented matrix product states"*, [arXiv:2412.17209](https://arxiv.org/abs/2412.17209) (2024).

## [PhasePoly (`:phasepoly`)](@id backend-phasepoly)

`PhasePoly` is an **optional package extension** powered by exact affine
basis permutations, phase polynomials, and quadratic Gauss sums. It adds zero
dependencies unless you install and load `PhasePoly`:

```julia
using Pkg
Pkg.add(url="https://github.com/aydindeger/PhasePoly.jl") # optional, once
using BlueTangle
import PhasePoly
```

Once `PhasePoly` is loaded, you can build a compact `PhasePolyState` with
`apply(ops; backend=:phasepoly, ...)` or `phasepoly(ops; N=..., max_h=..., max_cccz=..., max_ccs=...)`,
and pass `backend=:phasepoly` to `measure`, `sample`, and `shadow`:

```julia
using BlueTangle
import PhasePoly

# 1) Compact representation: exact native Pauli expectations & amplitudes
ops = [Op("H", 1), Op("H", 2), Op("CZ", 1, 2), Op("T", 1)]
compact = phasepoly(ops)
expect(compact, "X")                     # single-qubit expectations
correlation(compact, "X,Z", [1, 2])      # Gauss-sum Pauli expectation
amplitude(compact, [0, 1])               # single amplitude, exact global phase
@assert to_state(compact) ≈ apply(zero_state(2), ops)

# Single-branch Pauli queries scale to 80+ qubits without building a 2^N vector
large = phasepoly([Op("H", q) for q in 1:80])
@assert expect(large, "X", 71) ≈ 1
@assert correlation(large, "X,X", [2, 79]) ≈ 1

# 2) Unified backend=:phasepoly interface for sampling and shadows
m      = measure(ops, 2000; backend=:phasepoly, max_h=12, max_nonzeros=1 << 20)
draws  = sample(ops, 50; backend=:phasepoly)
rho_12 = shadow(ops, 2000; qubits=[1, 2], backend=:phasepoly)
```

Requesting `backend=:phasepoly` before loading `PhasePoly` raises an
`ArgumentError` reminding you to run `import PhasePoly`.

### Strengths and limitations (`:phasepoly`)

- **Strengths:** Tracks integer affine permutations and phase polynomials
  exactly, preserves **exact global phase** (in multiples of `π/16`), evaluates
  single-branch Pauli expectations and individual basis amplitudes at large `N`
  (`N = 80+`) without materializing `2^N` amplitudes, and treats leading and
  trailing Hadamard layers without branching.
- **Limitations:**
  - **Restricted gate set:** Supports only `H`, `X`, `Y`, `Z`, `S`, `SD`, `T`,
    `TD`, `HSP`, `HY`, `CX`/`CNOT`, `SWAP`, `CZ`, `CS`, `CSD`, `CT`, `CTD`,
    `CCZ`, `CCS`, `CCSdg`, single-qubit `P`/`U1`/`RZ` in integer multiples of
    `π/8`, and controlled phases (`CP`, controlled `P`/`U1`/`RZ`) in integer
    multiples of `π/4`. General continuous rotations (`RX(0.1)`, `RY`, `U3`),
    `CCX`, custom matrices, and variational functions raise `ArgumentError`.
  - **Noiseless & measurement-free only:** Gate noise (`Options(noise=...)`),
    readout noise (`Options(readout_noise=...)`), `OpQC`, mid-circuit
    measurements (`MX`/`MY`/`MZ`/`MR`), `RESET`, and `ifOp` all raise
    `ArgumentError`.
  - **Branching limits:** Interior (mid-circuit) `H` gates expand into `2^h`
    branches (guarded by `max_h = 12`). Level-4 (`√T`-rung) diagonal gates
    (`CT`, `CCS`, `CCCZ`) branch exponentially (`2^t`) over surviving `√T` terms
    in expectation queries (guarded by `max_cccz = 4`, `max_ccs = 8`).
  - **Sampling materializes the sparse vector:** `to_state`, `measure`, `sample`,
    `shadow`, and multi-branch expectations materialize a sparse statevector
    guarded by `max_nonzeros = 2^24` (and `N < 63`). Use `expect`, `correlation`,
    and `amplitude` directly on `phasepoly(ops)` to stay vector-free at large `N`.

## Choosing a backend

1. **Need exact ensemble averages of noise or unobserved measurements (`N ≲ 10–12`)?**
   Use `:density_matrix`.
2. **Need exact amplitudes, gradients, or reference checks (`N ≲ 18–20`)?**
   Use `:statevector`.
3. **Large `N` with 1D local / shallow entanglement and continuous rotation angles?**
   Use `:tensor` (`:mps`) and verify convergence in `maxdim` / `cutoff`.
4. **Large `N` with deep/nonlocal Clifford gates, stabilizer QEC, or Clifford+`T` circuits (`N_T ≲ N` or moderate magic)?**
   Use `:camps` (`CAMPSState`).
5. **Noiseless structured phase-polynomial / Clifford+`T` circuit with leading/trailing Hadamards and exact global phase at large `N`?**
   Install and `import PhasePoly`, then use `phasepoly(ops)` or `backend=:phasepoly`.
