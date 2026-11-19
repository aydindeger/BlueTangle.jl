# Quantum Error Correction

A stabilizer code stores logical information in the simultaneous +1 eigenspace
of commuting Pauli operators. An error can change their measurement outcomes
without revealing the encoded amplitudes. The resulting syndrome identifies a
candidate recovery operation.

BlueTangle constructs an encoding circuit from independent, commuting +1 Pauli
generators. Its built-in decoder searches for a matching **single-qubit Pauli**
error, which is sufficient for the example below. Correcting more general errors
requires a suitable decoder. The constructor stores the supplied distance `d`;
use `code_distance` to calculate distance from the stabilizers.

This section develops the workflow in three parts:

- **Overview and encoded algorithms** (this page): stabilizer codes, logical
  gates, and ideal correction rounds.
- [Noisy fault-tolerant gadgets](fault_tolerance.md): physical preparation,
  verified syndrome extraction, logical teleportation, logical T consumption,
  and readout.
- [Codes, distances, and magic-state factories](magic_factory.md): code
  comparisons, distance searches, noisy raw resources, and encoded distillation.

The [QEC example learning path](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/qec/README.md)
provides runnable steps in the order Bell pair → fault propagation → code
distance → noisy syndrome cycle → logical teleportation → distillation →
full noisy algorithm. The last example includes the encoded MPS resource
factory. Each step explains what to edit and how to interpret its output.

## Start developing encoded algorithms

For a concrete starting point, `steane_code(blocks)` builds independent
[[7,1,3]] Steane blocks: one logical qubit per seven consecutive physical
qubits. `steane_code(2)` therefore has 14 data qubits and two logical qubits.
The blocks form a product code with distance three for any number of blocks.
Concatenation would require encoding each physical qubit in a further code.

Logical `H`, `S`, `SD`, `X`, `Y`, `Z`, `CX`/`CNOT`, `CZ`, and `SWAP` have
transversal implementations. A one-qubit logical gate acts separately on the
physical wires; a two-qubit logical gate couples matching positions in distinct
blocks. Thus a fault at one physical gate cannot spread to several positions
within a block through these gates. The logical phase convention matters:
logical `S` uses physical `SD` on every wire, and logical `SD` uses physical `S`.
See [Gottesman's stabilizer-code thesis, §5.3](https://arxiv.org/html/quant-ph/9705052v1)
for the logical action of these transversal gates.

Use the existing `code.ops` interface to turn an ordinary logical circuit into
physical operations:

```@example qec_ft
using BlueTangle, LinearAlgebra

code = steane_code(2)
decoder = steane_decoder(code)
logical_ops = [Op("H", 1), Op("CX", 1, 2)]
physical_gates = code.ops(logical_ops)
expected = apply(zero_state(2), logical_ops)
encoded = apply(code.encode(), physical_gates)
@assert isapprox(abs2(dot(expected, code.decode(encoded))), 1; atol=1e-12)
(logical_qubits=code.k, physical_qubits=code.n, physical_gates=length(physical_gates))
```

This produces the logical Bell state `(|00⟩ + |11⟩)/sqrt(2)`. The circuit acts
on encoded data throughout; decoding at the end is a simulation diagnostic.
These mappings implement the Clifford part of an algorithm. Logical T/T† uses
`steane_t` with an explicit encoded magic resource and measurement feed-forward;
see [Noisy fault-tolerant gadgets](fault_tolerance.md) for this universal gate
set and physical Shor recovery. The [codes and factory guide](magic_factory.md)
adds Shor and Reed–Muller constructors, distance calculations, noisy raw
preparation, and the encoded 15-to-1 factory used by the full phase-estimation
example.

To simulate physical gate noise, use `apply(encoded, physical_gates;
noise=NoiseModel("depolarizing", 0.001))`. Noise follows the `noisy` flags on
the generated physical operations; setting `noisy=false` on a logical `Op`
does not change those flags.

### Repeated ideal correction rounds

`qec_round` applies ideal projectors directly in the data space without adding
ancillas. It returns a state containing only the original data qubits, so it
can be placed between logical gates and called repeatedly:

```@example qec_ft
# This is the output error produced by X_2 before the transversal CNOT.
damaged = apply(encoded, [Op("X", 2), Op("X", 9)])
first_round = qec_round(code, damaged; decoder=decoder)
second_round = qec_round(code, first_round.state; decoder=decoder)
@assert get_N(second_round.state) == code.n
@assert all(iszero, second_round.syndrome)
@assert isapprox(abs2(dot(expected, code.decode(second_round.state))), 1; atol=1e-12)
(first_syndrome=first_round.syndrome, second_syndrome=second_round.syndrome)
```

The Steane decoder treats X and Z components separately in each block. It
corrects one arbitrary single-qubit Pauli per block, including a Y error; a
fault propagated between two blocks can therefore be corrected even though the
whole register contains two errors. Supply `decoder=steane_decoder(code)`
explicitly: the generic default still searches for a single error across the
whole register. Higher-weight errors can produce the same syndrome and lead to
a logical error after recovery.

For statevectors and MPS, `qec_round` samples the ideal syndrome and returns
`(state=corrected_data, syndrome=bits)`. For density matrices it sums the
recovered sectors and returns `syndrome=nothing`; an ensemble has no single
record. It supports the same decoder callbacks as `code.correct`.
The `code.syndrome` and `code.correct(state)` paths retain their appended
syndrome ancillas, as explained below.

`qec_round` is an **ideal reference operation**: its stabilizer projectors and
recovery have no faults. The general synthesized encoder is also not a verified
fault-tolerant preparation circuit. The existing `code.ops_syndrome` uses a
bare ancilla for each check, and an ancilla fault can spread to multiple data
qubits. These distinctions matter when interpreting a successful simulation.

### Two runnable worked examples

From the repository root, run:

```sh
julia --project=. examples/qec/01_logical_bell.jl
julia --project=. examples/qec/02_fault_propagation.jl
```

- `01_logical_bell.jl` prepares an encoded Bell pair, injects a control-X error
  before logical CNOT, and shows its propagation to one qubit in each block.
  It corrects the errors, continues with phase and Hadamard gates, injects a
  second error, and corrects again. Both decoded fidelities should be 1 up to
  floating-point precision, with 14 data qubits after each round.
- `02_fault_propagation.jl` enumerates every nonidentity Pauli fault on the
  operands after each physical gate of logical `H₁, S₂, CX₁₂`. All 147 cases
  have trivial residual logical action after ideal recovery. It also exhibits
  an ancilla-X fault in bare-ancilla extraction that produces data `X₃X₄` and
  a residual logical X after recovery. This is a concrete failure of that
  extraction circuit's fault containment.

The second example uses `propagate_pauli(error, circuit; nqubits=code.n)` to
compute `U E U†`, ignoring its global phase. `pauli_syndrome(code, error)`
predicts the check outcomes without building a quantum state.
`logical_pauli(code, residual)` requires zero syndrome and returns binary
vectors `(x, z)` describing the remaining logical Pauli. Both vectors must be
zero for recovery to act trivially on arbitrary logical inputs. A zero
syndrome alone does not suffice. These helpers accept Pauli products as lists
of single-qubit `Op`s; unsupported gates in propagation are rejected.

The 147-case check covers one gate fault in that specified circuit with perfect
input preparation and final recovery. It omits idle, preparation, measurement,
and recovery faults and does not estimate a threshold. For physical extraction,
`shor_qec_round` now supplies verified ancillas and repeated full syndromes;
the [noisy-gadget guide](fault_tolerance.md) explains its fault model and tests.
Shor's [fault-tolerant error-correction construction, §4](https://arxiv.org/html/quant-ph/9605011v2)
explains why containing data errors and obtaining a reliable syndrome are
separate requirements. `steane_t` adds non-Clifford consumption, with magic-state
quality left as an explicit resource assumption.

## Encode an arbitrary logical qubit

The five-qubit code has four independent generators, leaving one logical qubit.
Commas separate the physical qubits in each generator.

```@example qec
using BlueTangle, LinearAlgebra

stabilizers = [
    "X,Z,Z,X,I",
    "I,X,Z,Z,X",
    "X,I,X,Z,Z",
    "Z,X,I,X,Z",
]
code = StabilizerCode(stabilizers; d=3)
logical = normalize(ComplexF64[1, 0.2 + 0.7im])
encoded = code.encode(logical)

@assert code.decode(encoded) ≈ logical
@assert all(s -> isapprox(expect(encoded, string_to_matrix(s)), 1), stabilizers)
(data_qubits=code.n, logical_qubits=code.k, syndrome_bits=code.m)
```

`code.encode()` prepares logical zero. Passing a statevector prepares its arbitrary
logical superposition. The physical logical-Pauli representatives are defined by
the constructed encoder; supplied logical-operation definitions take precedence.

## Inspect and correct a physical error

Apply a `Y` error to physical qubit 2, measure the generators, and inspect the
proposed correction:

```@example qec
damaged = apply(encoded, Op("Y", 2))
measured, bits, candidate = code.syndrome(damaged)
@assert candidate isa Op
(syndrome=bits, correction=(candidate.name, candidate.qubit),
  qubits_after_extraction=get_N(measured))
```

The bits follow the **input generator order**: 0 denotes eigenvalue +1 and 1
denotes eigenvalue −1. `code.syndrome` appends four ancillas and resets them after
reading the syndrome, so the returned state contains nine qubits.

The compact text diagram below shows just the injected error and the selected
Pauli correction. Syndrome extraction between them is implicit; encoding and
decoding are omitted from this small diagram.

```@example qec
plotq(QuantumOps[Op("Y", 2), candidate]; nqubits=5)
nothing # hide
```

```@example qec
corrected = code.correct(measured, bits)
recovered = code.decode(corrected)
fidelity_after_recovery = abs2(dot(logical, recovered))
@assert isapprox(fidelity_after_recovery, 1; atol=1e-12)
fidelity_after_recovery
```

The fidelity comparison ignores an irrelevant global phase. `decode` accepts
the five data qubits alone or the data plus reset syndrome ancillas. For a pure
state, it requires the discarded wires to be in zero before returning a decoded
logical state.

For a complete recovery cycle without inspecting the record:

```@example qec
recovered_directly = code.decode(code.correct(damaged))
@assert isapprox(abs2(dot(logical, recovered_directly)), 1; atol=1e-12)
```

A zero syndrome means no *detected* error. Some multi-qubit errors have a zero
syndrome or share a syndrome with a single-qubit error and can change the logical
state. A syndrome outside the decoder's search raises an error.

## Apply a logical gate

Here the qubit index refers to the logical register:

```@example qec
logical_X = code.apply(encoded, "X", 1)
decoded_X = code.decode(logical_X)
@assert decoded_X ≈ gate.X * logical
```

Use `apply(encoded, Op("X", physical_qubit))` when you intend a physical operation
instead. A physical error and a logical gate act on different objects even when
both are called `X`.

## Correct an ensemble of errors

A density matrix can describe a mixed logical input and a classical mixture of
physical errors. The following mixture contains no error, `X₁`, or `Y₄`, each of
which is correctable by this code:

```@example qec_density
using BlueTangle, LinearAlgebra
code = StabilizerCode([
    "X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"
]; d=3)
logical_rho = ComplexF64[0.6 0.2+0.1im; 0.2-0.1im 0.4]
encoded_rho = code.encode(logical_rho)
noisy_rho = 0.55encoded_rho +
            0.25apply(encoded_rho, Op("X", 1)) +
            0.20apply(encoded_rho, Op("Y", 4))
recovered_rho = code.decode(code.correct(noisy_rho))
@assert recovered_rho ≈ logical_rho
(trace=real(tr(recovered_rho)), error=norm(recovered_rho-logical_rho))
```

`code.correct(rho)` sums the projected syndrome sectors and their corrections
**deterministically**. Calling `code.syndrome(rho)` instead selects one Born
outcome and returns its conditional state and record:

```@example qec_density
conditional_rho, bits, candidate = code.syndrome(noisy_rho)
recovered_branch = code.decode(code.correct(conditional_rho, bits))
@assert recovered_branch ≈ logical_rho
@assert isapprox(real(tr(conditional_rho)), 1; atol=1e-12)
```

Density-matrix decoding traces out the other wires. This differs from pure-state
decoding, where discarding an entangled wire cannot produce another pure vector.

## Use an MPS without converting to a dense vector

Encoding, syndrome extraction, correction, and decoding also accept MPS inputs:

```@example qec_mps
using BlueTangle
code = StabilizerCode([
    "X,Z,Z,X,I", "I,X,Z,Z,X", "X,I,X,Z,Z", "Z,X,I,X,Z"
]; d=3)
logical_mps = plus_state(N_MPS(1))
encoded_mps = code.encode(logical_mps)
damaged_mps = apply(encoded_mps, Op("Y", 2))
recovered_mps = code.decode(code.correct(damaged_mps))
@assert recovered_mps isa it.MPS
@assert isapprox(abs2(inner(logical_mps, recovered_mps)), 1; atol=1e-10)
get_N(recovered_mps)
```

The MPS path samples a syndrome like a statevector trajectory. Its cost grows
with the bond dimensions generated by the encoding and recovery circuits, so
larger examples can be substantially more expensive.

## Reference

```@docs
StabilizerCode
steane_code
steane_decoder
qec_round
pauli_syndrome
propagate_pauli
logical_pauli
qec_state_prep
```
