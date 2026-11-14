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

Choose the output you need:

| Call | Returns |
|:--|:--|
| `apply(state, ops)` | Evolved state |
| `run(state, ops)` | Named tuple with `.state` and `.record` |
| `sample(state_or_ops, shots)` | Raw final outcomes as integers |
| `measure(state_or_ops, shots)` | Frequencies and expectation values |

Omit the shot count in `measure(state)` for exact probabilities, including on
MPS up to 16 qubits. Use `Measurement(draws, N)` to summarize existing raw draws.
Qubit indices start at 1, and qubit 1 is the most significant bit in statevectors
and integer measurement outcomes.

Circuit sampling starts each shot in the zero state. Gate noise is passed to
`apply` and `run(state, ops)` as `noise=...`, and to circuit sampling functions
as `options=Options(noise=...)`.
Operations run in their supplied order; no layout or SWAPs are added by default.
Use `route(ops, layout)` when connectivity matters. See the
[operation-list examples](examples.md#Operation-lists) for execution, sampling,
and saving results.

## Choose a workflow

- [Circuit examples](examples.md): gates, noise channels, measurements, plotting,
  and the distinction between an ensemble and a single trajectory. The
  [conditional measurement](@ref conditional-measurements) and
  [reset](@ref reset-and-records) sections explain outcome records and feedback.
  Later examples cover Hamiltonian evolution, entanglement,
  [partial trace](@ref partial-trace), Pauli twirling, and noise amplification.
- [Quantum Machine Learning](qml.md): VQE, held-out regression and XOR
  classification, gradient checks, noisy training, MPS truncation, and the
  complete QML API.
- [Quantum Error Correction](qec.md): stabilizer codes, logical gates, code
  distances, noisy syndrome extraction, and magic-state distillation. Worked
  examples progress from ideal correction to a complete noisy encoded algorithm.
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

PRIMA_jll version 0.7.1 can fail to load on Linux with glibc 2.41+
with `cannot enable executable stack`. See the
[upstream issue](https://github.com/libprima/PRIMA.jl/issues/46). Linux CI uses
Ubuntu 24.04 pending an upstream binary rebuild.

For a local checkout, start Julia in its root with `julia --project=.` and run
`using Pkg; Pkg.instantiate(); using BlueTangle`. Run the tests with
`julia --project=. -e 'using Pkg; Pkg.test()'`.

`plotq(ops)` prints a boxed text circuit without extra dependencies. For
graphical output, install and load `PyPlot`, then call
`plotq(ops; backend=:pyplot)` and `savefigure`. `plotqt` is also available for
text output. Loading `PhasePoly` enables `phasepoly` for its supported
structured noiseless circuits; see the
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

The default QEC decoder searches for a single physical Pauli error; the Steane
decoder treats each block independently. The QASM importer supports the
[documented OpenQASM 2 subset](qasm.md#Supported-subset). PhasePoly supports a
restricted gate family, and its cost can grow exponentially with interior
Hadamards or level-4 phase terms.
The relevant examples state these limits where they affect the calculation.

## Citation

Aydin Deger, “BlueTangle.jl: A Noisy and Dynamic Quantum Circuit Simulator”,
[repository](https://github.com/aydindeger/BlueTangle.jl/) (2024).
