# BlueTangle.jl

BlueTangle simulates quantum circuits with gates, noise, mid-circuit measurements,
and classical feed-forward. Statevectors, density matrices, and matrix product
states share the same circuit description. Stabilizer-code tools add encoding,
logical operations, syndrome extraction, and recovery.

## Start with a Bell pair

```@example home
using BlueTangle
ops = [Op("H", 1), Op("CX", 1, 2)]
state = apply(zero_state(get_N(ops)), ops)
@assert correlation(state, "Z,Z", [1, 2]) ≈ 1
@assert correlation(state, "X,X", [1, 2]) ≈ 1
(state=state, ZZ=correlation(state, "Z,Z", [1, 2]))
```

[`measure`](@ref) samples final outcomes; `run` returns the outcomes of explicit
mid-circuit measurements. Qubit indices start at 1, and qubit 1 is the most
significant bit in statevectors and integer measurement outcomes.

Use `apply(state, ops)` to evolve an existing state, or `measure(ops, shots)` to
prepare and sample repeated shots. Gate noise is passed to `apply` as
`noise=...`, and to sampling functions as `options=Options(noise=...)`.
Operations run in their supplied order; no layout or SWAPs are added by default.
Use `route(ops, layout)` when connectivity matters. See the
[operation-list migration notes](examples.md#Operation-lists-and-migration) for
replacements for the removed wrapper API.

## Choose a workflow

- [Circuit examples](examples.md): gates, noise channels, measurements, plotting,
  and the distinction between an ensemble and a single trajectory. The
  [conditional measurement](@ref conditional-measurements) and
  [reset](@ref reset-and-records) sections explain outcome records and feedback.
  Later examples cover Hamiltonian evolution, entanglement,
  [partial trace](@ref partial-trace), Pauli twirling, and noise amplification.
- [Quantum Machine Learning](qml.md): VQE, held-out regression and XOR
  classification, gradient crosschecks, noisy training, MPS truncation, and the
  complete QML API.
- [Quantum error correction](qec.md): encode an arbitrary logical qubit in the
  five-qubit code, inspect a syndrome, correct a physical error, and decode;
  repeat with density matrices or an MPS.
- [OpenQASM import and export](qasm.md): read strings or files, use whole-register
  operations, preserve idle qubits, and diagnose unsupported input.
- [Function reference](functions.md): signatures and detailed API documentation.

## Installation

Julia 1.10 or later is required. Install from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/aydindeger/BlueTangle.jl")
using BlueTangle
```

The current PRIMA_jll 0.7.1 dependency can fail to load on Linux with glibc 2.41+
with `cannot enable executable stack`. See the
[upstream issue](https://github.com/libprima/PRIMA.jl/issues/46). Linux CI uses
Ubuntu 24.04 pending an upstream binary rebuild.

For a local checkout, start Julia in its root with `julia --project=.` and run
`using Pkg; Pkg.instantiate(); using BlueTangle`. Run the tests with
`julia --project=. -e 'using Pkg; Pkg.test()'`.

`plotq(ops)` prints a boxed text circuit without extra dependencies. For
graphical output, install and load `PyPlot`, then call
`plotq(ops; backend=:pyplot)` and `savefigure`. `plotqt` is also available for text. Loading `PhasePoly`
enables `phasepoly` for its supported structured noiseless circuits; see the
[optional backend example](examples.md#Optional-PhasePoly-backend).

Install optional packages in the environment where you use them:

```julia
using Pkg
Pkg.add("PyPlot")
Pkg.add(url="https://github.com/aydindeger/PhasePoly.jl")
```

## Scope

The package is experimental. Density-matrix evolution sums noise and unobserved
measurement branches, whereas vectors and MPS states represent individual
trajectories. MPS truncation, finite-shot sampling, and Trotter discretization
introduce different approximations; choose and check their settings separately.

The QEC decoder searches single-Pauli errors. The QASM importer supports an
explicit OpenQASM 2 subset. PhasePoly supports a restricted gate family, and its
cost can grow exponentially with interior Hadamards or level-4 phase terms.
The relevant examples state these limits where they affect the calculation.

## Citation

Aydin Deger, “BlueTangle.jl: A Noisy and Dynamic Quantum Circuit Simulator”,
[repository](https://github.com/aydindeger/BlueTangle.jl/) (2024).
