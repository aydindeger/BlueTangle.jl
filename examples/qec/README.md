# Learning fault-tolerant algorithms with BlueTangle

These examples build up from an encoded Bell pair to a complete noisy
Clifford+T computation with its own magic-state factory. Start with the short
examples to understand which errors are corrected, which are detected and
discarded, and which can survive as logical errors. Run the full factory after
those distinctions are familiar.

Run commands from the repository root. Each script prints its results to the
terminal and checks its central claims with assertions. They perform local
simulations and submit no hardware jobs. Julia's first run also compiles the
package, so it takes longer than later runs.

## Suggested order

The numbering preserves the existing examples; the learning order is
**01 → 02 → 04 → 05 → 06 → 03**.

| Example | Question it answers | Simulation and cost |
|:--|:--|:--|
| [01: Logical Bell pair](01_logical_bell.jl) | Can a propagated physical fault be corrected without losing logical entanglement? | 14-qubit data statevector; ideal encoding and correction |
| [02: Fault propagation](02_fault_propagation.jl) | Why are transversal gates safe against one gate fault, while a bare syndrome ancilla can be dangerous? | Fast binary Pauli propagation; no large statevector |
| [04: Codes and distances](04_codes_and_distance.jl) | What distinguishes Shor, Steane, and Reed–Muller codes, and how is distance verified? | Bounded support search and statevectors up to 15 qubits |
| [05: Noisy syndrome cycle](05_noisy_syndrome_cycle.jl) | How do verified cats and repeated syndromes protect a noisy recovery cycle? | One Steane block plus five Shor ancillas; small statevectors |
| [06: Distillation walkthrough](06_distillation_walkthrough.jl) | Why does 15-to-1 reject one bad raw input but sometimes accept three? | Six-wire physical reference circuit and exact probability formulas |
| [03: Full factory and phase estimation](03_noisy_phase_estimation.jl) | How does an algorithm obtain and consume a newly distilled resource on every trial? | Up to 42 entangled MPS wires; 47 physical wires including factory ancillas; more expensive |

The first five are small standalone calculations. Example 03 simulates many
physical preparations, gates, checks, retries, and readouts for each logical T.
Its cost depends on noise, rejected attempts, and MPS entanglement. The physical
wire count and simulation representation are different quantities: independent
preparation circuits are simulated separately until their registers interact.

## 1. Encoded Bell pair

```sh
julia --project=. examples/qec/01_logical_bell.jl
```

Two Steane blocks store two logical qubits. An X on physical wire 2 propagates
through a transversal CNOT to `X₂X₉`. There are now two physical errors in the
whole register, but only one in each block. The blockwise decoder can correct
both. The printed Bell fidelity and the fidelity after another round of gates
and correction should be one, up to floating-point error.

Extend `next_gates` to try further Clifford operations after recovery. The
fault-insertion section explicitly executes H followed by CNOT, so changing
`logical_circuit` also requires updating that schedule and its predictions.
Change `fault` and the later Z error to explore propagation; update assertions
that explicitly predict `X₂X₉` when changing the worked fault. Encoding,
`qec_round`, and the final decoder are ideal here. Example 05 introduces noisy
syndrome extraction.

## 2. Fault propagation

```sh
julia --project=. examples/qec/02_fault_propagation.jl
```

The script checks all 147 single gate-location Pauli faults in its short
transversal Clifford circuit. After ideal blockwise recovery, it finds zero
logical failures for that finite set. It then puts one X on a bare syndrome
ancilla after its second data interaction. That error spreads to two remaining
data wires, and the decoder adds a third X: the residual is a logical X even
though its syndrome is zero.

Change `logical_circuit` to audit another Clifford circuit, and update the
expected fault count if the number of physical gates changes. This example
makes the distinction between **zero syndrome** and **correct logical state**
concrete. It also motivates the verified-cat circuit in example 05. Its audit
does not include faulty preparation, readout, or recovery.

## 3. Code comparison and distance

```sh
julia --project=. examples/qec/04_codes_and_distance.jl
```

The searched distances are Shor `dX=dZ=3`, Steane `dX=dZ=3`, and Reed–Muller
`dX=7,dZ=3`. Thus all three have `d=3`. The output also lists a
minimum-weight logical Pauli. A separate search stopped after weight two
returns `distance=nothing` and `lower_bound=3`, establishing `d≥3` while leaving
the exact distance undetermined. An arbitrary complex logical input is then
recovered after one physical Y, and the Reed–Muller transversal logical T is
checked.

Change `max_weight` or `max_candidates` in `code_distance` to control search
cost. Inspect `exact`, `lower_bound`, `searched_through`, and `truncated` before
interpreting a result. `code.d` stores the value supplied to the constructor;
`code_distance` calculates distance from the stabilizer checks. The search
cost grows exponentially with code size in general.

Shor's **nine-qubit code** is distinct from Shor's **verified-cat extraction
method**. The current noisy extraction applies the latter to Steane blocks;
it does not convert them to the nine-qubit code.

## 4. Noisy syndrome extraction

```sh
julia --project=. examples/qec/05_noisy_syndrome_cycle.jl
```

This example first physically prepares logical plus, corrects Y on data wire
3, and prints the syndrome history. It then deliberately flips one cat
readout: the first record disagrees with the second, so a third matching
record is needed before recovery. Finally, a noisy logical-zero preparation,
transversal H, noisy correction, and decoded X readout form a complete
Clifford cycle. A changed syndrome record calls for another round;
`QECProtocolFailure` means the bounded protocol did not finish and is not a
logical measurement outcome.

The Y₃ record is `[1,0,1,1,0,1]` in this code's check order, repeated twice.
The deliberate readout fault instead gives `[1,0,0,0,0,0]`, then two all-zero
records. A syndrome bit of one means a check measured eigenvalue −1; these
records locate faults without directly measuring the logical amplitudes.

Edit `gate_error`, `prep_error`, and `measurement_error` near the top, each
initially `5e-4`. The script prints the records needed to follow the cycle.
The detailed derivation and callback convention
are in the [noisy-gadget guide](../../docs/src/fault_tolerance.md).

## 5. Magic-state distillation

```sh
julia --project=. examples/qec/06_distillation_walkthrough.jl
```

The six-wire reference prepares and consumes fifteen physical raw magic
states. A single deliberately bad raw state makes the first attempt fail;
its checks are `[0,0,0,1]`. The next clean attempt succeeds, so the count is
two attempts and 30 raw states. Defects on columns 1, 2, and 3 instead pass
because their binary labels XOR to zero. The accepted output has a logical Z
error: some faulty states can pass all four checks.

The final probability table compares the leading `35p³` term with the exact
conditional output error, acceptance probability, and mean raw-state cost
`15/P_accept`. These formulas assume independent raw Z errors and perfect
Clifford operations, with unlimited independent retries for the mean cost.
Vary the raw error probability to inspect that limit. To study Clifford faults,
pass a gate-noise model to `distill15`. This unencoded reference is for learning
and verification. Example 03 uses
`steane_distill15` to protect the distillation circuit itself.

## 6. Phase estimation with a resource factory

```sh
julia --project=. examples/qec/03_noisy_phase_estimation.jl
```

Each trial prepares a fresh noisy raw-resource supply, runs encoded 15-to-1
distillation until acceptance or its retry bound, then prepares computational
`|+⟩`. The accepted resource implements T; noisy Shor recovery and decoded X
or Y readout finish the trial. There is no ideal encoded magic input or ideal
recovery in this operational path.

For a perfect T, the readout means are `cos(π/4)` and `sin(π/4)`. The script
estimates those means and then the phase. The default two shots per basis show
the workflow and produce estimates with large uncertainty. Increase `shots`
for statistics; each additional shot needs its own complete factory. Probability estimates
use completed trials, while resource counts include rejected attempts and
aborted trials. Approximate 95% Wilson intervals describe sampling uncertainty
of the zero probabilities and remain broad even when both shots agree. They
do not measure physical-model or MPS truncation error. The attempt records
show which four-bit checks passed or caused a retry.

The main settings at the top are:

| Setting | Meaning |
|:--|:--|
| `shots` | Independent full trials per X or Y basis |
| `gate_error` | Depolarizing probability per operand after a physical gate |
| `prep_error` | X probability after a fresh physical zero preparation |
| `measurement_error` | Probability of a flipped physical measurement bit |
| `raw_error` | Additional Z probability after each raw physical T |
| `factory_max_attempts` | Maximum complete factory attempts per trial |
| `cutoff` | MPS singular-value truncation; numerical accuracy, not physical noise |

The factory runs before the computational data is prepared, giving a physical
peak of 47 wires. Keeping that extra Steane block alive during distillation
would raise the count to 54. Idle noise is not inserted implicitly. Reduce
`cutoff` and check convergence before interpreting quantitative MPS results.

The [factory guide](../../docs/src/magic_factory.md) derives the protocol and
resource counts and describes the single-fault checks and multi-fault
simulations. These examples do not establish a threshold or a noisy advantage
over an unencoded algorithm.
