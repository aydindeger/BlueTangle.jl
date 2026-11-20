# Clifford-Augmented MPS (CAMPS)

Matrix product states (`:tensor` / `:mps`) compress states with low spatial
entanglement, but even pure Clifford circuits (and stabilizer error-correcting
codes) can generate volume-law entanglement whose MPS bond dimension grows
exponentially with system size.

**Clifford-Augmented Matrix Product States (CAMPS)** separate stabilizer
entanglement from non-Clifford ("magic") entanglement by writing an `N`-qubit
pure state as

```math
|\Psi\rangle = C |\phi\rangle,
```

where:

- `C` is an `N`-qubit Clifford unitary tracked in `O(N²)` bits by a
  [`CliffordTableau`](@ref), and
- `|\phi\rangle` is a disentangled non-Clifford `ITensorMPS.MPS` core.

In `CAMPSState`, every Clifford gate updates `C` in `O(N)` time without
touching `|\phi\rangle` or increasing its bond dimension. Non-Clifford gates
pull their Pauli generators through `C` (`P̃ = C† P C`), collapse them onto
untouched `|0⟩` sites when available, apply the rotation to `|\phi\rangle`, and
sweep two-qubit Clifford disentanglers across adjacent bonds to move newly
created stabilizer entanglement back into `C`.

## Quick start

Use [`camps`](@ref) to evolve `|0…0⟩` directly into a [`CAMPSState`](@ref), or
pass `backend=:camps` (alias `:clifford_mps`) to `measure`, `sample`, `run`, or
`shadow`:

```@example camps_intro
using BlueTangle

ops = [
    Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3), Op("CX", 3, 4),
    Op("T", 2), Op("RY(0.6)", 3), Op("CZ", 1, 4),
]

st = camps(ops; N=4)
sv = apply(zero_state(4), ops)

@assert maxlinkdim(st) == 1
@assert fidelity(st, sv) ≈ 1
@assert isapprox(expect(st, "Z"), expect(sv, "Z"); atol=1e-12)
@assert isapprox(correlation(st, "X,Z", [1, 3]), correlation(sv, "X,Z", [1, 3]); atol=1e-12)
(bond_dim=maxlinkdim(st), fidelity=fidelity(st, sv), Z_expect=round.(expect(st, "Z"); digits=4))
```

You can also initialize `CAMPSState(N)` explicitly, wrap an existing `ITensorMPS.MPS`
with `CAMPSState(psi)`, or use the same `apply`, `run`, and `*` syntax as
statevectors and MPS:

```@example camps_intro
st0 = CAMPSState(4)
st1 = apply(st0, ops)          # non-mutating: st0 remains |0000⟩
st2 = Op("S", 1) * st1
@assert maxlinkdim(st0) == 1
@assert expect(st0, "Z") ≈ ones(4)
@assert fidelity(st2, apply(sv, Op("S", 1))) ≈ 1
```

## How CAMPS works

Every `CAMPSState` holds two fields: `state.tableau::CliffordTableau` and
`state.mps::it.MPS`. The tableau stores the inverse-conjugated Pauli basis
`px[j] = C† Xⱼ C` and `pz[j] = C† Zⱼ C` for `j = 1:N`.

1. **Clifford gates (`O(N)` cost, zero MPS work):**
   Applying a Clifford unitary `U` updates `C ← U C`, which transforms each
   tracked row via `(UC)† P_j (UC) = C† (U† P_j U) C`. The MPS core `|\phi\rangle`
   and its bond dimension are completely untouched. This also applies to
   parameterized gates at exact multiples of `π/2` (such as `RX(pi/2)`,
   `RZZ(pi/2)`, `RXX(pi/4)`, or `CP(pi)`), and to Pauli noise channels.

2. **Pauli pullback for non-Clifford rotations:**
   A Pauli rotation `exp(-i (θ/2) P)` acting on `C |\phi\rangle` equals
   `C exp(-i (θ/2) P̃) |\phi\rangle`, where `P̃ = C† P C` is computed directly
   from the tableau in `O(N)` bitvector operations.

3. **Free-zero Pauli reduction (exact, zero bond growth):**
   Before applying `exp(-i (θ/2) P̃)` to `|\phi\rangle`, CAMPS identifies any
   untouched product `|0⟩` sites in `|\phi\rangle` (where `⟨φ|Z_q|φ⟩ ≈ 1`).
   - If `P̃` has `X` or `Y` on a free zero `q₀`, controlled Cliffords with
     control `q₀` leave `|\phi\rangle` identically unchanged while conjugating
     `P̃` down to a weight-1 operator `±X_{q₀}`. The rotation then acts locally
     on site `q₀` with zero bond-dimension increase.
   - If `P̃` has only `I` or `Z` on free zeros, those `Z` factors act as `+1`
     and are dropped.

   Similarly, when a projective Pauli measurement collapses onto a free-zero
   pivot `q₀`, the outcome is sampled in `O(N²)` tableau operations and a
   single-qubit Clifford restores site `q₀` to `|0⟩`, recovering exact
   Gottesman–Knill evolution at `maxlinkdim == 1` for Clifford circuits with
   mid-circuit measurements and resets.

4. **Two-qubit Clifford disentangling sweeps:**
   When a multi-qubit `P̃` must be applied via an MPS direct sum
   `cos(θ/2)|\phi\rangle - i sin(θ/2) P̃ |\phi\rangle`, CAMPS sweeps across
   adjacent bonds `(b, b+1)` testing the 19 left cosets of
   `C₁^{⊗2} \ C₂` (identity, 9 CNOT-class cosets, and 9 double-CNOT/iSWAP-class
   cosets). Whenever a two-qubit Clifford `D` lowers the bond von Neumann
   entropy, CAMPS updates `|\phi\rangle ← D |\phi\rangle` and `C ← C D†`,
   keeping `C |\phi\rangle` invariant while shrinking the MPS core.

## Bond-dimension compression vs. standard MPS

In nonlocal or deep Clifford circuits doped with non-Clifford (`T`) gates, a
plain MPS hits the volume-law ceiling `χ = 2^{N/2}`, whereas CAMPS absorbs the
Clifford entanglement into `state.tableau`. Here, a 10-qubit nonlocal Clifford
circuit doped with 10 `T` gates reaches `χ = 32` in plain MPS (and in CAMPS
with `disentangle=false`), while default CAMPS stays at `χ = 1` with machine
precision fidelity:

```@example camps_compress
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

sv = apply(zero_state(N), ops)
mps_plain = apply(zero_state(N_MPS(N)), ops; cutoff=1e-12)
st_camps = camps(ops; N=N, cutoff=1e-12, disentangle=true, nsweeps=2)
st_nodis = camps(ops; N=N, cutoff=1e-12, disentangle=false)

@assert fidelity(st_camps, sv) ≈ 1
@assert maxlinkdim(st_camps) < maxlinkdim(mps_plain)
(
    plain_mps_chi = maxlinkdim(mps_plain),
    camps_nodis_chi = maxlinkdim(st_nodis),
    camps_chi = maxlinkdim(st_camps),
    fidelity = fidelity(st_camps, sv),
)
```

### Truncation and disentangling controls

All `CAMPSState` evolution entry points (`apply`, `run`, `camps`, `measure`,
`sample`, `shadow`, and `qec_round`) accept the following keyword arguments:

| Keyword | Default | Role |
|:--|:--|:--|
| `cutoff` | `1e-12` | Singular-value truncation threshold for MPS direct sums and sweeps. |
| `maxdim` | `nothing` | Optional hard cap on the core MPS bond dimension (`nothing` = uncapped). When set, intermediate direct sums retain up to `2 * maxdim` before disentangling compresses back down to `maxdim`. |
| `disentangle` | `true` | Enable free-zero Pauli reduction and two-qubit Clifford disentangling sweeps. Set `false` to disable disentangling (useful for benchmarking). |
| `nsweeps` | `2` | Maximum number of left/right two-qubit disentangling sweeps per multi-qubit non-Clifford operation or projection. |

Once the number of independent non-Clifford gates exceeds `N`, magic
entanglement begins to grow in `state.mps`. When running with a finite `maxdim`
cap, verify convergence by comparing fidelities or observables across budgets,
just as for `:tensor`:

```@example camps_compress
# Add extra Clifford+T layers beyond the N=10 free-zero budget
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
@assert maxlinkdim(st_chi2) <= 2
@assert maxlinkdim(st_chi8) <= 8
(
    chi2_fid = round(fidelity(st_chi2, sv_more); digits=5),
    chi8_fid = round(fidelity(st_chi8, sv_more); digits=5),
)
```

## Supported operations

CAMPS supports BlueTangle's full circuit model:

- **Zero-MPS-cost Clifford gates:** `I`, `X`, `Y`, `Z`, `H`, `S`, `SD`, `SX`
  (`XSQRT`), `SXDG` (`SXD`, `XSQRTD`, `XSQRTDG`), `HSP`, `HY`, `H2`, `CI`, `CX`
  (`CNOT`), `CZ`, `CY`, `SWAP`, `ISWAP`, `FSWAP`, `ECR`, plus controlled Paulis
  and any parameterized gate at an exact multiple of `π/2`.
- **Single-qubit non-Clifford & arbitrary `U(2)` gates:** `T`, `TD`, `P`, `U1`,
  `U2`, `U3`, `RX`, `RY`, `RZ`, and custom 2×2 unitary matrices (pure Z, X, and
  Y rotations use a single Pauli rotation; general unitaries use
  `zyz_decomposition`).
- **Two-qubit non-Clifford & arbitrary `U(4)` gates:** `CS`, `CSD`, `CT`,
  `CTD`, `CH`, `CP`, `RXX`, `RYY`, `RZZ`, `RXY`, `GIVENS`, `FSIM`, `SYC`,
  `SWAPA`, and arbitrary custom 4×4 unitaries via exact Kraus–Cirac magic-basis
  KAK decomposition.
- **Controlled and three-qubit gates:** arbitrary controlled single-qubit
  unitaries (`Op(U, q; control=c)` via Barenco ZYZ decomposition), plus `CCZ`,
  `CCX`, `CCY`, and `CSWAP`.
- **Variational gates:** function-backed `Op(name, f, q)` evaluated via `pars=...`.
- **Mid-circuit measurements, reset, and feedback:** `MX`, `MY`, `MZ`, `MR`,
  `born_measure_Z`, `Op("RES", q)` / `Op("RESET", q)`, and conditional `ifOp`
  (including measurements inside branches).
- **Noise and channels:** `OpQC` and `NoiseModel` (Pauli channels update the
  tableau directly with zero MPS cost; coherent rotation channels use the
  unitary path; non-unitary channels such as `amplitude_damping` and
  `phase_damping` compute exact Born weights from `K†K` and apply the sampled
  Kraus operator in the pulled-back Pauli basis).
- **Composite blocks (`OpF`):** operation-list blocks, state callbacks returning
  a `CAMPSState`, and 1- or 2-qubit dense matrices.

## Circuit sampling, records, noise, and classical shadows

Use `backend=:camps` (or `backend=:clifford_mps`) with `measure`, `sample`,
`run`, and `shadow`. All `Options` fields—gate noise, readout noise, readout
mitigation, measurement basis, Pauli twirling, `layout` routing, and
`density_matrix=true` (for `N <= 16`)—work identically to `:statevector` and
`:tensor`:

```@example camps_api
using BlueTangle, LinearAlgebra

ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("CX", 2, 3)]

# Exact expectations from a CAMPS state, and finite-shot circuit sampling
m = measure(ops, 2000; backend=:camps, cutoff=1e-12, maxdim=16)
draws = sample(ops, 10; backend=:camps)

# Mid-circuit measurement, reset, and conditional feedback
dynamic_ops = [
    Op("H", 1),
    Op("MZ", 1; if0=Op("I", 2), if1=Op("X", 2)),
    Op("RES", 1),
    Op("MZ", 1),
]
traj = run(CAMPSState(2), dynamic_ops)
@assert traj.record[2] == 0 # qubit 1 was reset to |0⟩

# Noisy trajectories compared against the exact density matrix
nm = NoiseModel(Noise1("amplitude_damping", 0.25), Noise2("depolarizing", 0.1))
rho_exact = apply(to_rho(zero_state(3)), ops; noise=nm)
m_noisy = measure(ops, 3000; backend=:camps, options=Options(noise=nm))
@assert isapprox(m_noisy.expect, expect(rho_exact, "Z"); atol=0.08)

# Subsystem classical shadow without materializing a full statevector
rho_13 = shadow(ops, 3000; qubits=[1, 3], backend=:camps, cutoff=1e-12)
@assert size(rho_13) == (4, 4)
(noisy_Z = round.(m_noisy.expect; digits=3), exact_Z = round.(expect(rho_exact, "Z"); digits=3))
```

When the MPS core is the product state `|0…0⟩` (for example after any Clifford
circuit, stabilizer QEC round, or Clifford circuit with Pauli noise),
`sample_bit` and `sample` automatically take a pure bitvector stabilizer-sampling
path with zero ITensor contractions per shot. Outcomes wider than a machine
integer use `BigInt`.

## Stabilizer QEC at bond dimension 1

Standard MPS simulation of multi-block stabilizer codes incurs substantial bond
dimension from the encoding circuits and transversal entangling gates alone.
Because stabilizer encoders, transversal Clifford gates, Pauli errors, and
stabilizer syndrome projections all live in the Clifford group, `CAMPSState`
executes them in `tableau` without growing `mps`:

```@example camps_qec
using BlueTangle

# 3 Steane [[7,1,3]] blocks = 21 physical qubits
code = steane_code(3)
decoder = steane_decoder(code)

# Encode |000⟩_L and prepare a logical GHZ state using transversal Cliffords
enc = code.encode(CAMPSState(3))
enc = code.apply(enc, "H", 1)
enc = code.apply(enc, "CX", 1, 2)
enc = code.apply(enc, "CX", 2, 3)
@assert get_N(enc) == 21
@assert maxlinkdim(enc) == 1

# Inject one physical Pauli error into each of the three blocks
damaged = apply(enc, [Op("X", 3), Op("Y", 11), Op("Z", 19)])

# Ideal stabilizer projection and block-wise decoding stays at chi = 1
recovered = qec_round(code, damaged; decoder=decoder)
@assert maxlinkdim(recovered.state) == 1
@assert fidelity(recovered.state, enc) ≈ 1

# Decode back to a 3-qubit logical CAMPSState
logical = code.decode(recovered.state)
@assert get_N(logical) == 3
@assert correlation(logical, "X,X,X", [1, 2, 3]) ≈ 1
@assert correlation(logical, "Z,Z", [1, 2]) ≈ 1
(
    physical_qubits = get_N(enc),
    bond_dim_after_qec = maxlinkdim(recovered.state),
    syndrome = recovered.syndrome,
    logical_XXX = correlation(logical, "X,X,X", [1, 2, 3]),
)
```

Non-Clifford logical states (such as arbitrary single-qubit `U3` states or
`T`-doped logical circuits) also encode and decode cleanly with `code.encode`,
`code.apply`, `code.syndrome`, `code.correct`, `qec_round`, and `code.decode`;
the Steane encoder itself still adds zero MPS bond dimension.

## Inspecting the Clifford frame and converting states

Use [`clifford_ops`](@ref) to synthesize an equivalent `O(N²)` sequence of
elementary Clifford gates for `state.tableau`, and [`to_state`](@ref),
[`to_rho`](@ref), or [`to_MPS`](@ref) to materialize the physical state when
needed:

```@example camps_convert
using BlueTangle

ops = [Op("H", 1), Op("T", 1), Op("CX", 1, 2), Op("S", 2)]
st = camps(ops; N=2)
sv = apply(zero_state(2), ops)

c_ops = clifford_ops(st)
@assert all(op -> !op.noisy, c_ops)
@assert same_up_to_global_phase(to_state(st), sv)[1]
@assert to_rho(st) ≈ sv * sv'
@assert abs2(inner(to_MPS(st; cutoff=1e-12), sv)) ≈ 1
(synthesized_clifford_gates = length(c_ops), mps_core_chi = maxlinkdim(st))
```

!!! note "Global phase semantics"
    `CliffordTableau` tracks the Clifford frame `C` as an element of the
    projective Clifford group (up to a single global phase `e^{iα}`).
    Consequently, `to_state(state)` and `to_MPS(state)` equal the exact
    physical state up to a single global phase (`same_up_to_global_phase`),
    while `to_rho(state)`, `expect`, `correlation`, `fidelity`, `sample`,
    `sample_bit`, `sample_exact`, `measure`, `prob`, `run`, `shadow`, and
    `qec_round` are exact and independent of global phase.

## Choosing a backend and limitations

| Backend | Best suited for | Limitations |
|:--|:--|:--|
| `:statevector` | Exact pure-state trajectories up to `~16–20` qubits | `O(2^N)` memory and time |
| `:density_matrix` | Exact ensemble averages of noise & measurements up to `~10–12` qubits | `O(4^N)` memory |
| `:tensor` (`:mps`) | Large `N` with low spatial entanglement (1D area-law circuits) | Bond dimension grows exponentially under deep/nonlocal Clifford entanglement |
| `:camps` | Large `N` with deep Clifford entanglement / QEC and low-to-moderate non-Clifford (`T`) count | Magic-dominated circuits (`N_T ≫ N`) still grow the core MPS bond dimension exponentially if uncapped |
| `phasepoly` | Noiseless structured Clifford+`T`/phase-polynomial circuits | Optional extension; restricted gate set; exponential in interior Hadamards |

Keep the following operational limits in mind when using `:camps`:

- **Magic-dominated circuits:** CAMPS removes stabilizer entanglement, not
  non-stabilizer ("magic") entanglement. Once the number of independent
  non-Clifford gates substantially exceeds `N`, the core MPS bond dimension
  still grows unless capped with `maxdim` (which introduces truncation error).
- **Exact enumeration limit:** `sample_exact(state)` and `measure(state)`
  (without a shot count) materialize `2^N` amplitudes and therefore default to
  `max_qubits = 16`. Pass a shot count (`measure(state, shots)`) for larger
  registers.
- **Re-absorbing the Clifford frame:** `to_MPS(state)` applies
  `clifford_ops(state)` back onto `state.mps`. If `C` is deeply entangling, the
  resulting plain MPS can have a much larger bond dimension than `state.mps`;
  prefer evaluating `expect`, `correlation`, `sample`, and `qec_round` directly
  on `CAMPSState`.
- **Full-register dense matrices:** Dense `OpF(name, mat)` matrices above 2
  qubits are unsupported on `CAMPSState`; decompose them into 1- and 2-qubit
  operations.

See also the runnable walkthrough script
[`examples/camps_clifford_t.jl`](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/camps_clifford_t.jl).

## [API reference](@id CAMPS-API-reference)

```@docs
CAMPSState
CliffordTableau
camps
clifford_ops
```

`CAMPSState` also integrates with BlueTangle's shared `apply`, `run`,
`measure`, `sample`, `shadow`, `expect`, `correlation`, `fidelity`,
`born_measure_Z`, `to_state`, `to_rho`, `to_MPS` (documented in
[Functions](functions.md)), and `StabilizerCode` / `qec_round` (documented in
[Quantum Error Correction](qec.md)).
