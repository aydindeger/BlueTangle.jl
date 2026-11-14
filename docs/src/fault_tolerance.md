# Noisy fault-tolerant gadgets

BlueTangle supplies a concrete distance-three Steane workflow: physical logical
zero/plus preparation, verified Shor syndrome extraction, transversal Clifford
gates, logical T/T† consumption, and destructive logical X/Z readout. These are
small statevector or MPS simulations of physical circuits. They use the same
`steane_code` convention as the [QEC tutorial](qec.md).

`steane_t` consumes an explicitly supplied encoded magic state. The small
consumption check below uses an ideal input to isolate that gadget; the full
phase-estimation script obtains every resource from `steane_distill15`, which
simulates noisy physical preparation and encoded 15-to-1 distillation. A
logical resource error can still become a logical computational error. The
[codes and factory guide](magic_factory.md) explains how the factory filters
such errors, its resource costs, and the assumptions behind its protection.
For the runnable sequence from ideal code behavior to the complete noisy
algorithm, follow the [QEC learning path](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/qec/README.md).

## Why the syndrome circuit is different

A single ancilla controlling four data interactions can spread one X fault to
several data qubits. The Shor construction distributes those interactions over
four ancillas in the cat state `(0000 + 1111)/sqrt(2)`. Each cat qubit interacts
with exactly one data qubit, and the parity of their X measurements gives the
stabilizer eigenvalue.

The cat itself must be checked: an error while making it can create correlated
X errors. Here its fixed preparation order is `H₁, CX₁₂, CX₁₃, CX₁₄`. Before
touching the data, a fifth ancilla measures `Z₁Z₂`. Outcome one discards that
attempt and starts a fresh cat. For this specific order, the check detects every
dangerous cat X pattern produced by one preparation fault, up to the harmless
`XXXX` cat stabilizer. Changing the preparation order requires a new audit.

A phase error in a cat or an incorrect readout can still report a wrong
syndrome. Therefore, `shor_qec_round` measures **complete syndromes until two
consecutive records agree**, then applies the recovery. It does not take an
independent majority vote for each check. This is the distance-three instance
of the traditional repeated-syndrome Shor protocol described by
[Tansuwannont, Pato and Brown](https://arxiv.org/abs/2208.05601).

```@example noisy_qec
using BlueTangle, LinearAlgebra
code = steane_code()
logical = normalize(ComplexF64[1, 0.2+0.7im])
damaged = apply(code.encode(logical), Op("Y", 3))
result = shor_qec_round(code, damaged)
@assert get_N(result.state) == 7
@assert abs2(dot(logical, code.decode(result.state))) ≈ 1
(syndromes=result.history, cat_attempts=result.attempts)
```

This zero-noise check uses the actual ancilla circuit. With noise enabled, a
fault-tolerant correction round may leave one correctable physical error.
That is consistent with the distance-three requirement; exact return to the
code space is not guaranteed after a fault near the end of the round.

## Noise, records, and failure bounds

```julia
noise = NoiseModel("depolarizing", 0.0005)
result = shor_qec_round(code, damaged;
    noise=noise, prep_error=0.0005, measurement_error=0.0005)
```

The parameters refer to distinct physical locations:

| Parameter | Meaning |
|:--|:--|
| `noise` | The existing gate-noise model, after cat, verification, data-coupling, readout-basis, and physical recovery gates |
| `prep_error` | Probability of X immediately after each fresh physical `\|0⟩` preparation |
| `measurement_error` | Probability of a flipped physical Z outcome, represented by X just before measurement |
| `fault` | Optional callback injecting a specified Pauli on a location's operands |

For `NoiseModel("depolarizing", p)`, a two-qubit gate applies independent Pauli
noise with probability `p` on **each operand**. Its probability of any error is
`1-(1-p)^2`. Exhaustive fault tests also include every correlated two-operand
Pauli at one gate. Idle errors, leakage, faults in classical processing, and
correlations between distinct locations are not included implicitly.

`result.history` contains complete syndrome records; `result.attempts[r][j]`
counts cat attempts for check `j` in round `r`. `result.trace` records physical
locations and measurement outcomes. The trace is local to each function call.
The default limits are six full rounds and twenty cat attempts per check.
Exhaustion raises `QECProtocolFailure`; simulations must count such aborted
shots separately. A failed preparation attempt never interacts with the data.

`shor_syndrome_circuit(code, check)` exposes the exact per-check schedule.
For reproducible fault tests, a callback receives `index`, `stage`, `round`,
`check`, `attempt`, `kind`, and `op`. It returns one Pauli `Op` or a vector of
them on that location's operands. Gate/preparation faults occur after the
operation; measurement faults occur before it. For example:

```@example noisy_qec
# This creates a dangerous cat error before verification, so the cat is rejected.
fault = loc -> loc.stage == :cat_encode && loc.round == 1 && loc.check == 1 &&
    loc.attempt == 1 && loc.op.name == "CX" && loc.op.target_qubit == code.n+3 ?
    Op("X", code.n+1) : Op[]
checked = shor_qec_round(code, code.encode(logical); fault=fault)
@assert checked.attempts[1][1] == 2
@assert abs2(dot(logical, code.decode(checked.state))) ≈ 1
checked.attempts[1]
```

## Physical preparation and logical readout

`steane_prepare(:zero)` starts with seven physical zero states and uses the
verified Shor circuit to project onto the code space. Logical Z is initially
positive; random X-check sectors are corrected with Z operations. `:plus`
adds physical H gates before extraction. This path includes preparation faults
and avoids the general synthesized encoder.

```@example noisy_qec
prepared = steane_prepare(:plus)
@assert abs2(dot(plus_state(1), code.decode(prepared.state))) ≈ 1
readout = steane_measure(code, prepared.state; basis=:X)
@assert readout.outcome == 0
@assert length(readout.state) == 1  # the only block has been consumed
(physical_bits=readout.bits, corrected_bits=readout.corrected_bits)
```

Readout measures all seven physical wires and classically decodes their Hamming
syndrome. The parity of the corrected word is the logical outcome. One physical
input error **or** one measurement error is tolerated. For logical X, physical
H gates precede Z measurements. There is no inverse encoder in this path.
`steane_readout(bits)` exposes the classical decoder separately. With several
blocks, `steane_measure(...; block=b)` removes only that block.

`steane_prepare(...; backend=:tensor)` creates an MPS. The Shor, injection, and
readout functions also accept MPS inputs directly. Tensor simulation cost can
still grow with entanglement; these examples do not establish scalability.
`steane_t` and `steane_measure` accept an optional MPS `cutoff`; the encoded
factory sets it explicitly. This numerical truncation is separate from the
physical noise model and should be varied when checking simulation convergence.

## Worked cycle: correction, disagreement, and readout

Before running the resource factory, follow one small syndrome cycle:

```sh
julia --project=. examples/qec/05_noisy_syndrome_cycle.jl
```

The script first prepares logical plus through the physical preparation
circuit, inserts Y on data wire 3, and corrects it with actual Shor extraction.
The syndrome histories identify the physical defect, while the decoded
fidelity checks preservation of the logical state. Decoding here diagnoses a
simulation; it is not a step in the noisy operational circuit.
For Y₃ the six-bit syndrome is `[1,0,1,1,0,1]` in this constructor's check
order, and two identical records precede recovery. A one denotes a measured
check eigenvalue of −1.

Next, one deliberately flipped cat readout changes the first complete
syndrome. The next fault-free record disagrees with it, so recovery must wait
for another matching record. The printed history makes the stopping rule
visible: one misleading record is insufficient to trigger recovery.

Finally, a fresh noisy logical-zero preparation, transversal H, noisy Shor
round, and decoded logical X measurement form a complete Clifford cycle.
`gate_error`, `prep_error`, and `measurement_error` are editable near the top;
all default to `5e-4`. A readout bit is one stochastic trajectory, while
`QECProtocolFailure` is a separately reported heralded abort. The script does
not infer an error rate or threshold from that one trajectory. The
[learning path](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/qec/README.md)
places this cycle between code-distance checks and magic-state distillation.

## T and T† by consuming magic states

Let `|A⟩ = (|0⟩ + exp(iπ/4)|1⟩)/sqrt(2)`. For data `α|0⟩ + β|1⟩`, a CNOT
from data to `|A⟩`, followed by Z measurement of the resource, gives one of
two data branches:

```math
m=0:\quad T|\psi\rangle/\sqrt{2},\qquad
m=1:\quad e^{i\pi/4}T^\dagger|\psi\rangle/\sqrt{2}.
```

Applying S for outcome one makes both branches T, up to an irrelevant global
phase. The encoded gadget uses seven transversal CNOTs, fault-tolerant logical
Z readout of the resource, and transversal conditional S. With `dagger=true`,
supply `|A†⟩` and the conditional correction is S†. This is the resource-state
approach to non-Clifford computation; resource purification is a separate task
from consumption, as in [Bravyi and Kitaev](https://arxiv.org/abs/quant-ph/0403025).
The following zero-noise check isolates consumption; `code.encode(target)`
here is a diagnostic input, not the preparation used by the full algorithm.

```@example noisy_qec
target = apply(plus_state(1), Op("T", 1))
magic = code.encode(target)  # supplied ideal resource; NOT a physical factory
prepared = steane_prepare(:plus)
injected = steane_t(code, prepared.state, magic)
corrected = shor_qec_round(code, injected.state)
@assert abs2(dot(target, code.decode(corrected.state))) ≈ 1
(injection_outcome=injected.outcome, remaining_data_qubits=get_N(corrected.state))
```

The resource is a required argument. There is no default that silently prepares
it using the unverified encoder. `steane_t` accepts physical gate noise,
measurement noise, and deterministic fault callbacks; no ideal correction is
hidden inside it. Separate gadgets with `shor_qec_round` to prevent physical
errors from accumulating. A logical Z error on the resource remains a logical
error after stabilizer correction, and is covered by a negative-control test.

## Worked algorithm: estimating a π/4 phase

Run from the repository root:

```sh
julia --project=. examples/qec/03_noisy_phase_estimation.jl
```

Each trial of this small Ramsey/phase-estimation algorithm first prepares and
distills a fresh encoded magic state with the noisy factory. It then prepares
logical `|+⟩`, consumes the accepted resource to apply T, corrects the data,
and reads logical X or Y. The ideal probabilities are

```math
P_X(0)=\frac{1+\cos\phi}{2},\qquad
P_Y(0)=\frac{1+\sin\phi}{2},\qquad \phi=\pi/4.
```

The sample averages estimate the two quadratures, and `atan(sin_phi, cos_phi)`
estimates the phase. The default two shots per basis demonstrate the complete
workflow, including resource accounting, rather than providing precise phase
statistics. The script reports completed and aborted shots, approximate 95%
Wilson intervals for the completed-shot zero probabilities, factory attempts,
and raw-state consumption. Those intervals remain broad even if both default
shots agree; they describe binomial sampling uncertainty, not physical-model
or MPS truncation error. Factory
counts include rejected attempts and work in aborted trials; probability and
phase estimates are conditioned on completed shots. Edit `shots`, the physical
error probabilities, `raw_error`, the retry bound, and `cutoff` at the top.

The factory runs before the computational data is prepared. Its physical peak
is 47 wires, including five Shor ancillas, while the entangled MPS contains at
most 42 wires. The script uses no ideal encoded magic input or ideal recovery.
`raw_error` adds an independent Z after each raw physical T; gate, preparation,
and readout errors are simulated separately throughout the factory and
computation. The printed `distillation15_statistics` values assume perfect
Clifford operations and serve only as a reference. The finite noisy run is not
a threshold estimate or evidence of an advantage over an unencoded circuit.

## What the tests establish

An independent binary propagation test follows the actual Shor schedule and
adaptive stopping rule. It checks all 2,160 single-location Pauli faults on a
clean encoded input, all 21 incoming single-qubit Paulis with perfect recovery,
and all possible initial projection sectors of zero/plus preparation. It also
compares representative histories with real statevector trajectories. The
single-fault cases require at most four syndrome rounds.

These recovery claims use the default Steane decoder. A custom decoder is
checked for syndrome consistency, which alone cannot certify its logical action.

The injection tests compare T/T† with the intended logical action, including an
entangled input, both feed-forward branches, physical CNOT/input/resource/
readout/correction faults, and a logical-resource-error negative control.
These checks establish the stated finite distance-three fault model; they do
not extend it to arbitrary multi-fault trajectories or unmodelled noise.
The [factory guide](magic_factory.md) distinguishes this inner-code protection
from the outer 15-to-1 filter, and derives its ideal-Clifford error formula.

```@docs
shor_qec_round
shor_syndrome_circuit
steane_prepare
steane_t
steane_measure
steane_readout
QECProtocolFailure
```
