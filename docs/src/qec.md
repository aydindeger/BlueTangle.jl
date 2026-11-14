# Quantum error correction

A stabilizer code stores logical information in the simultaneous +1 eigenspace
of commuting Pauli operators. An error can change their measurement outcomes
without revealing the encoded amplitudes. The resulting syndrome identifies a
candidate recovery operation.

BlueTangle constructs an encoding circuit from independent, commuting +1 Pauli
generators. Its built-in decoder searches for a matching **single-qubit Pauli**
error. This is sufficient for the example below; it is not a decoder for arbitrary
multi-qubit errors, and setting `d=3` does not independently certify a code's distance.

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
state, it checks that discarded wires are in zero rather than silently throwing
away unresolved information.

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

The MPS path samples a syndrome like a statevector trajectory. It can still
become expensive when the encoding and recovery circuits generate large bond
dimensions; a small example succeeding is not a scaling guarantee.

## Reference

```@docs
StabilizerCode
qec_state_prep
```
