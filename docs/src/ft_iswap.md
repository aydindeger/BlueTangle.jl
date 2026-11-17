# A fault-tolerant logical iSWAP

This page builds a small *early fault-tolerant* gadget: a logical iSWAP between two
error-detecting code blocks, made only from the gates that superconducting processors
with √iSWAP couplers execute natively. A runnable version with assertions is
[`examples/qec/07_fault_tolerant_iswap.jl`](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/qec/07_fault_tolerant_iswap.jl).

## What "fault tolerant" means here

A code of distance 2 cannot correct errors, but it detects any single-qubit error.
A circuit is **fault tolerant** at this level if every *single fault* is either detected
or harmless. A single fault is one wrong Pauli after any gate, a wrong preparation, or a
flipped measurement. Detected shots are discarded.
Then an undetected logical error needs at least two faults, so its probability scales as
`p²` for a physical error rate `p`, instead of `p` for unprotected qubits.

The obstacle is that one faulty two-qubit gate can damage *both* qubits it acts on.
Inside a single block this can create an undetectable logical error from one fault.
Transversal gates avoid this: each physical gate connects one qubit of block A with the
corresponding qubit of block B, so a fault leaves at most one error per block.

## The code: [[6,2,2]] on a row of six qubits

[`c6_code`](@ref) is Knill's C6 code. Its checks act on the same two supports for X and Z:

| position | 1 | 2 | 3 | 4 | 5 | 6 |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| check a (X and Z) | ● | | ● | ● | ● | |
| check b (X and Z) | | ● | ● | | ● | ● |
| logical qubit 1 (X_L = XXX, Z_L = ZZZ) | ● | ● | ● | | | |
| logical qubit 2 (X_L = XXX, Z_L = ZZZ) | | | | ● | ● | ● |

The logical supports have **odd weight**, so each block behaves like a pair of
Steane-like qubits:

* H on every wire is logical H on both qubits;
* SD on every wire is logical S on both qubits (S on every wire is logical S†);
* ISWAP on the six pairs of corresponding wires of two blocks is logical iSWAP† on both
  logical pairs, and with Z on the first block before and after it is logical iSWAP.

The popular [[4,2,2]] code cannot do this. All of its logical operators have even X and Z
weight. Its transversal gates therefore preserve the parity of the number of Y's in every
logical Pauli, while iSWAP maps X₁ to Z₁Y₂. A physical iSWAP between two data qubits of one
[[4,2,2]] block is a logical iSWAP, but a single XX, YY or ZZ fault after it is an
undetectable logical error. [[6,2,2]] is the smallest code of this type with odd-weight logicals.

## Building and checking the circuit

```julia
using BlueTangle

# block 1 in |-,+>_L (X basis), block 2 in |0,1>_L (Z basis), then two logical iSWAPs
circuit = vcat(c6_prepare(:X; bits=(1, 0), block=1), c6_prepare(:Z; bits=(0, 1), block=2),
               c6_iswap(), c6_iswap(), c6_readout(:X; block=1), c6_readout(:Z; block=2))
psi = apply(zero_state(12), circuit)
shot = c6_decode(int2bin(sample(psi, 1)[1], 12), (:X, :Z))
shot.accepted   # true: every parity check holds
shot.logical    # [0, 1, 0, 1] (A1, A2, B1, B2): iSWAP² = Z⊗Z flips the two X-basis values
```

[`c6_prepare`](@ref) uses H on two wires and eight nearest-neighbour CNOTs in four layers.
An exhaustive search over short CNOT circuits on a row shows that seven CNOTs suffice to
prepare logical |00⟩, and that eight are needed if no single fault may leave an undetectable
logical error. [`c6_single_faults`](@ref) checks the whole circuit. It propagates every
single fault with [`propagate_pauli`](@ref), which now supports `ISWAP`, and classifies the result:

```julia
audit = c6_single_faults(circuit)
audit.harmful    # 0: no single fault is both undetected and harmful
```

## Noise: p versus p²

[`sample_pauli_frames`](@ref) samples depolarizing faults through a Clifford circuit
and returns the final Pauli frame of every shot. This is exact for Pauli noise and much
faster than a state-vector simulation. Example 07 uses it to compare eight logical iSWAPs,
with every gate noisy at rate `p` and ideal preparation and readout:

| p | bare qubits | one [[4,2,2]] block | two [[6,2,2]] blocks (kept shots) |
|:--|:--|:--|:--|
| 1e-3 | 8.0e-3 | 1.6e-3 | 1.5e-5 |
| 4e-3 | 3.1e-2 | 6.6e-3 | 2.2e-4 |
| 1.6e-2 | 1.2e-1 | 2.7e-2 | 3.9e-3 |
| log-log slope | 0.98 | 1.03 | 2.03 |

The bare and [[4,2,2]] errors grow linearly in `p`. The transversal [[6,2,2]] gadget's
undetected error grows quadratically. The price is the discarded shots, and the circuits
have no mid-circuit checks. Two separately detectable errors made at different times can
therefore combine into an undetected one, and the kept error grows faster than linearly
in the number of gates.

## Reference

```@docs
c6_code
c6_logicals
c6_prepare
c6_iswap
c6_readout
c6_decode
c6_single_faults
sample_pauli_frames
```
