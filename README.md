# BlueTangle.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deger.xyz/BlueTangle.jl)
[![Build Status](https://github.com/aydindeger/BlueTangle.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aydindeger/BlueTangle.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Documentation

[**BlueTangle.jl Documentation**](https://deger.xyz/BlueTangle.jl)

## Introduction

BlueTangle is a Julia toolkit for quantum circuits with unitary gates, noise,
mid-circuit measurements, and classical feed-forward. The same operations work
with statevectors, density matrices, and matrix product states. It also includes
stabilizer-code encoding and recovery, Hamiltonian evolution, variational
circuits, and an optional PhasePoly backend for structured noiseless circuits.

Start with the examples below, then see the
[QEC tutorial](https://deger.xyz/BlueTangle.jl/qec/) and
[OpenQASM guide](https://deger.xyz/BlueTangle.jl/qasm/).
The package is experimental; the supported semantics and limits are documented
alongside each workflow.

BlueTangle is developed by [Aydin Deger](https://deger.xyz).

## Installation
Installing BlueTangle is straightforward. You can install it through the Julia REPL using the following command:

```julia
using Pkg
Pkg.add(url="https://github.com/aydindeger/BlueTangle.jl")
```

Alternatively, type `] add https://github.com/aydindeger/BlueTangle.jl` directly in the Julia REPL.

To load the package in your Julia environment, use:

```julia
using BlueTangle
```

On Linux with glibc 2.41 or newer, the current PRIMA_jll 0.7.1 dependency can
prevent loading with `cannot enable executable stack`. This is tracked in an
[open upstream issue](https://github.com/libprima/PRIMA.jl/issues/46); Linux CI
uses Ubuntu 24.04 while awaiting an upstream binary rebuild.

For a local checkout, start Julia from the package root with `julia --project=.` and run:

```julia
using Pkg
Pkg.instantiate()
using BlueTangle
```

To run the tests from the package root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Documentation dependencies are kept in `docs/Project.toml`. To build the documentation locally:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.

## Basic usage

For more examples and advanced usage, please refer to the [**BlueTangle.jl Documentation**](https://deger.xyz/BlueTangle.jl)

### Creating and Analyzing a GHZ Circuit

This example demonstrates the creation of a Greenberger–Horne–Zeilinger (GHZ) state circuit, including the steps for measurement and calculating correlations.

```julia
# Quantum operations to prepare the GHZ state
hadamard = Op("H", 1)
cnot1 = Op("CNOT", 1, 2)
cnot2 = Op("CNOT", 2, 3)
ops = [hadamard, cnot1, cnot2]

# Evolve a state directly, or prepare fresh shots from the operation list.
state = apply(zero_state(3), ops)
shots = 1000
measurement = measure(ops, shots)

# Output measurement details
println("Expectation values:", measurement.expect)
println("Total magnetization moments:", measurement.mag_moments)

# Calculate and print correlations (e.g., ⟨Z₁Z₂⟩)
correlations = correlation(measurement, [1, 2])
println("Correlations:", correlations)
```

## Backends and measurement

Julia 1.10 or later is required. The same operation list works across backends:

```julia
ops = [Op("H", 1), Op("CX", 1, 2)]
measure(ops, 1000; backend=:statevector)
measure(ops, 1000; backend=:tensor, cutoff=1e-12)
measure(ops, 1000; backend=:density_matrix)
```

`measure` samples final outcomes. `run` returns the explicit mid-circuit measurement
records. Applying operations to a density matrix sums measurement and noise branches;
`apply(...; track_measurements=true)` instead samples measurement outcomes and
returns the conditional density matrix with its record.

The [measurement and reset examples](docs/src/examples.md) explain Pauli outcome
bits, conditional branches, records from measurements inside a branch, and MPS
truncation. `Op("RESET", q)` resets without recording its internal outcome;
`ifOp("MZ", q, "I", "X")` with tracking retains the pre-reset bit. The examples
also explain why explicit feedback gates can receive gate noise while standalone
reset is ideal. See [reduced states and partial trace](docs/src/examples.md) for
the distinction between keeping qubits in the vector API and discarding
subsystems in the density-matrix API.

Readout mitigation uses signed inverse-response weights, which can be negative
at finite shot counts.

MPS `measure(psi)` defaults to 1000 shots. Use `measure(psi, shots)` to set the
sample count, or `measure(psi, -1)` for exact probabilities on up to 16 qubits.
`sample_exact(psi; max_qubits=...)` explicitly overrides that enumeration limit.
Sampled results above 16 qubits carry an empty `density_matrix` placeholder.
MPS conversion and integer outcomes follow BlueTangle's convention: qubit 1 is
the most significant bit. Outcomes wider than a machine integer use `BigInt`.

Classical shadows can reconstruct a selected subsystem directly with MPS:

```julia
rho12 = shadow(ops, 1000; qubits=[1, 2], backend=:mps,
               cutoff=1e-12, maxdim=128) # A 4×4 estimate, without a full statevector.
```

Omitting `qubits` requests the full matrix; the default backend is `:statevector`.
See the [MPS shadow example](docs/src/examples.md#two-qubit-mps-estimate-versus-the-exact-state)
for an exact two-qubit reference and error versus shot count. Finite-shot shadow
estimates need not be positive; MPS truncation is a separate source of error.

Parameterized gates accept arithmetic and numeric functions, such as
`Op("RY(sin(pi/4))", 1)`. Gate strings do not execute Julia code. For custom
parameterized operations, pass a function directly: `Op("rotation", theta ->
_RX(theta), 1)`.

## Operation lists and optional routing

`apply(state, ops)` evolves the supplied state in list order. The state type
selects the backend; pass `noise=...` for gate noise and
`track_measurements=true` to return `(state, record)`.

`measure(ops, shots; options=Options(noise=noise))` prepares independent noisy
shots. Sampling one noisy state repeatedly would instead sample one realized
trajectory. Density-matrix evolution sums noise branches and can reuse the
resulting ensemble. Readout noise, mitigation, measurement bases, and per-shot
Pauli twirling are configured through `Options`.

Routing is explicit: `routed = route(ops, layout)` inserts the needed SWAPs and
returns another operation list. `measure`, `run`, and `shadow` also accept an
optional `layout` keyword. Without one, the given gates are applied directly.
Use `get_stats(ops)`, `get_layers(ops)`, and `plotq(ops)` for inspection.

Save operations and settings with
`save_data(ops, measurement; options=options, path="results.ser")`;
`load_data` returns their `ops`, `options`, `measurement`, and `other` fields.

## Optional packages

`plotq(ops)` prints a boxed circuit diagram without a plotting dependency.
`plotqt` provides direct text rendering. For graphical output,
install and load PyPlot, then select it explicitly:

```julia
using Pkg
Pkg.add("PyPlot") # once, in the environment used for plotting
using BlueTangle, PyPlot
fig = plotq(ops; backend=:pyplot)
savefigure(fig, "circuit.png")
```

All documentation figures can be regenerated from the repository root with
`julia --project=/path/to/plot-env docs/generate_figures.jl`, using an environment
containing BlueTangle and PyPlot. The script includes the finite-size MPS
trajectory experiment and its standard-error bars.

For supported noiseless circuits, install PhasePoly separately and load it to
enable the optional adapter:

```julia
using Pkg
Pkg.add(url="https://github.com/aydindeger/PhasePoly.jl") # optional, once
using BlueTangle
import PhasePoly
state = phasepoly([Op("H", 1), Op("H", 2), Op("CZ", 1, 2), Op("T", 1)])
expect(state, "X")                     # native compact Pauli calculation
correlation(state, "X,Z", [1, 2])
to_state(state)                        # explicitly materialize a sparse vector
```

The adapter preserves BlueTangle's qubit ordering and global phase. Unsupported
gates and noise are rejected explicitly. PhasePoly's interior-Hadamard and
level-4 phase branching can still be exponential; it is not a general replacement
for the other backends. See `?phasepoly` for the supported gate set and limits.

## Error correction and evolution

A five-qubit code can preserve an arbitrary logical state through a single
physical Pauli error:

```julia
using BlueTangle, LinearAlgebra
code = StabilizerCode([
    "X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"
]; d=3)
logical = normalize(ComplexF64[1, 0.2 + 0.7im])
encoded = code.encode(logical)
damaged = apply(encoded, Op("Y", 2))
recovered = code.decode(code.correct(damaged))
@assert abs2(dot(logical, recovered)) ≈ 1
```

`code.syndrome` exposes the syndrome record and proposed correction. Vector and
MPS recovery samples a syndrome; `code.correct(rho)` on a density matrix sums
syndrome sectors deterministically. The decoder searches single-Pauli errors;
it is not a general decoder for arbitrary multi-qubit errors. See the
[QEC tutorial](https://deger.xyz/BlueTangle.jl/qec/) for logical gates, mixed states,
and native MPS examples.

For encoded algorithm experiments, `steane_code(2)` supplies two Steane blocks
with transversal logical Clifford gates. Compile a logical circuit with
`code.ops(...)`, interleave gates with data-only `qec_round` calls, and use
`steane_decoder(code)` to recover one Pauli error per block. The worked examples
are runnable from the repository root:

```sh
julia --project=. examples/qec/01_logical_bell.jl
julia --project=. examples/qec/02_fault_propagation.jl
```

The first creates an encoded Bell state and corrects a propagated fault. The
second audits every gate-local Pauli fault in a short transversal circuit and
demonstrates a harmful fault in the bare-ancilla syndrome circuit. Recovery in
`qec_round` is ideal; these are building blocks for fault-tolerance studies,
not a complete fault-tolerant preparation, recovery, or universal gate set.

`hamiltonian_exp` approximates `exp(-im*time*H)`, matching `hamiltonian` and the
exact `evolve(state, H, time)` helper. It uses a shorter final Trotter step to reach
the requested time, including when `time` is not an integer multiple of `dt`.

## Quantum Machine Learning

Define a scalar loss and train its parameter vector directly:

```julia
ops = [Op("RY", _RY, 1)]
energy(parameters) = real(expect(variational_apply(parameters, ops), gate.Z))
result = train(energy, fill(0.3, nparameters(ops)); iterations=30)
result.loss, result.parameters, result.history
```

`train` returns the same fields for each optimizer. `variational_apply` accepts
an `init` state, so one circuit can evaluate many feature inputs. The
[QML section](https://deger.xyz/BlueTangle.jl/qml/) brings together VQE, regression
and XOR classification with held-out checks, gradient validation, deterministic
noise, MPS truncation, and the QML API reference.
`AnsatzOptions`, `VQA`, and `VQE` remain available for existing code.

## Import and export OpenQASM

```julia
source = """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0],q[1];
"""
ops = from_qasm(source; is_file=false)
state = apply(zero_state(get_N(ops)), ops)
to_qasm(ops, "bell.qasm")
restored = from_qasm("bell.qasm"; is_file=true)
@assert apply(zero_state(get_N(restored)), restored) ≈ state
```

The importer supports an OpenQASM 2 gate subset, indexed and whole-register
operands, multiple quantum registers, reset, and measurements. Errors include
the source line. Custom gate definitions and classical conditions must be
expanded before import; OpenQASM 3 is not supported. See the
[OpenQASM guide](https://deger.xyz/BlueTangle.jl/qasm/) for measurement records,
safe arithmetic, and interoperability limits.

## Constructing operations

`Op` is the common entry point for local gates, conditional measurements,
channels, and complete-register blocks:

```julia
Op("H", 1)
Op("MZ", 1; if0="I", if1="X")
Op(Noise1("bit_flip", 0.1), 1)
Op("bell block", [Op("H", 1), Op("CX", 1, 2)])
```

The returned objects retain distinct semantics: a gate, measurement branch,
Kraus channel, or complete-register transformation. Existing `ifOp`, `OpQC`, and
`OpF` constructors remain available. A parameterized *local* gate takes a qubit
index, for example `Op("rotation", theta -> _RX(theta), 1)`; `Op("transform", f)`
is a function of the whole state.

Noise amplification is explicit: copy an operation list, call
`cnot_amplifier!`, and measure the resulting circuit. Simulation never adds
noise-scaling circuits automatically. The [manual folding and extrapolation example](docs/src/examples.md#13-manual-noise-amplification-and-extrapolation)
states the assumptions behind that workflow.

## Citation 

For citation in research work, please use the following reference:

```
Aydin Deger, "BlueTangle.jl: A Noisy and Dynamic Quantum Circuit Simulator", https://github.com/aydindeger/BlueTangle.jl/ (2024).
```

