# Codes, distances, and magic-state factories

A fault-tolerant computation needs two kinds of protection. Error correction
keeps an already encoded state usable after physical faults. Magic-state
distillation checks several imperfect non-Clifford resources and keeps an
output only when those checks agree with the desired state. An undetected
logical resource error survives ordinary stabilizer correction, while faults
in the distillation circuit require protected Clifford operations.

BlueTangle supplies three small codes, a bounded distance calculation,
physical noisy raw-state preparation, a six-wire reference distillation circuit,
and an encoded 15-to-1 factory. The encoded factory produces a resource that
can be consumed by `steane_t`, as shown below. The reference circuit is
useful for understanding the protocol and checking formulas; its Clifford gates
are unprotected.

The [QEC learning path](https://github.com/aydindeger/BlueTangle.jl/blob/main/examples/qec/README.md)
puts runnable code-distance and six-wire distillation walkthroughs before the
more expensive encoded factory. Start there to reproduce the mechanisms one
at a time.

## Shor, Steane, and Reed–Muller codes

The notation `[[n,k,d]]` means `n` physical qubits, `k` encoded qubits, and
distance `d`: the smallest number of physical positions on which a nontrivial
logical Pauli can act. For a CSS code, `dX` and `dZ` restrict the search to
X-only and Z-only operators, and `d = min(dX,dZ)`.

| Constructor | Parameters | `dX`, `dZ` | Role in these examples |
|:--|:--|:--|:--|
| `shor_code()` | `[[9,1,3]]` | 3, 3 | A concrete repetition-based error-correcting code |
| `steane_code()` | `[[7,1,3]]` | 3, 3 | The blocks used for noisy Clifford gates, extraction, and readout |
| `reed_muller_code()` | `[[15,1,3]]` | 7, 3 | The code underlying 15-to-1 magic-state distillation |

Shor's nine-qubit code first groups the wires into three triples. Within each
triple, two Z-pair checks locate a bit flip. Two X checks, each spanning a pair
of triples, compare their phases. A Z error on any member of one triple has the
same effect on the encoded information; their differences are stabilizers.
For example, `Z₁Z₂` has weight two but acts trivially in the code space, so its
existence does not make the distance two.

Shor syndrome extraction is a verified-cat measurement method.
`shor_qec_round(steane_code(), state)` uses it to measure the stabilizers of
seven-qubit Steane data. The current noisy extraction implementation is
specialized to Steane blocks; `qec_round` provides ideal correction for all
three code constructors.

The Reed–Muller constructor labels its fifteen physical columns by the nonzero
four-bit words, in increasing integer order. The four coordinate rows give its
X checks; the ten Z checks span the binary kernel of those rows together with
the all-ones row. Their independent ranks give `15-4-10=1` logical qubit.
Physical `T` on all fifteen wires implements logical T† in this convention.
Accordingly, `rm.ops("T", 1)` uses fifteen physical TD gates.

## Code distance

The constructor stores the declared distance in `code.d`. The function
`code_distance` calculates distance by searching Pauli supports in increasing
weight. A candidate must commute with
every stabilizer and anticommute with at least one logical X/Z representative.
The second condition excludes stabilizers such as Shor's `Z₁Z₂`.

```@example code_distance
using BlueTangle
rm = reed_muller_code()
result = code_distance(rm; max_weight=7)
@assert result.exact && result.distance == 3
@assert (result.dx, result.dz) == (7, 3)
@assert all(iszero, pauli_syndrome(rm, result.witness))
(distance=result.distance, dX=result.dx, dZ=result.dz,
 witness=[(op.name, op.qubit) for op in result.witness])
```

The Z witness is on columns 1, 2, and 3. Their labels satisfy
`0001 XOR 0010 XOR 0011 = 0000`, so that Z product commutes with all four X
checks. Its odd weight anticommutes with logical X on all fifteen wires, hence
it is a logical operator. Exhausting weights one and two proves that its weight
three is minimal.

A bounded search may prove only a lower bound:

```@example code_distance
limited = code_distance(rm; max_weight=2)
@assert limited.distance === nothing
@assert limited.lower_bound == 3
(exact=limited.exact, lower_bound=limited.lower_bound,
 searched_through=limited.searched_through, stop_reason=limited.stop_reason)
```

This result establishes `d ≥ 3`. Equality additionally requires a weight-three
logical operator. `max_candidates` also limits the work, and `truncated=true`
identifies a budget interruption.
Only fully completed weights contribute to a lower bound. CSS codes need only
pure X and pure Z searches; a general stabilizer code uses X/Y/Z assignments on
each support. The method is bounded enumeration, so a large search can still be
expensive.

Run the executable code comparison, including ideal error correction and a
Reed–Muller transversal T check, from the repository root:

```sh
julia --project=. examples/qec/04_codes_and_distance.jl
```

## Noisy magic-state preparation

The resource consumed by a logical T gate is

```math
|A\rangle = T|+\rangle
          = \frac{|0\rangle+e^{i\pi/4}|1\rangle}{\sqrt{2}}.
```

`raw_magic_state` prepares a physical zero, applies H and T, and samples the
requested preparation and gate noise. An additional `phase_error=p` applies Z
after T with probability `p`; it lets us isolate the standard distillation
input model from other gate errors.

```@example factory_reference
using BlueTangle, LinearAlgebra
target = apply(plus_state(1), Op("T", 1))
raw = raw_magic_state()
@assert abs2(dot(target, raw.state)) ≈ 1
[(event.stage, event.op.name) for event in raw.trace]
```

The encoded raw-preparation path applies the physical Steane encoder and Shor
recovery. Faults in that encoder can produce logical errors that survive
recovery. These remaining logical resource errors enter the distillation
protocol.

## The 15-to-1 distillation circuit

The reference circuit uses five persistent wires and one reusable raw-resource
wire. Persistent wire 1 is the output; wires 2–5 hold four checks. For column
`j` with binary label `(b₁,b₂,b₃,b₄)`, define the Z product

```math
P_j = Z_1 Z_2^{b_1} Z_3^{b_2} Z_4^{b_3} Z_5^{b_4},\qquad j=1,\ldots,15.
```

Prepare the five wires in plus states. For each column, parity CNOTs compute
that Z product onto wire 1, one raw magic state implements a T on that parity,
and the CNOTs are undone. The net operation is
`exp(-iπ Pⱼ/8)`, up to a global phase. The circuit consumes fifteen raw
resources, each with its own preparation and noise.
This is the compressed Reed–Muller circuit described in
[Litinski, §3.1](https://arxiv.org/abs/1808.02892).

For example, column 3 is `0011`, so its support is `[1,4,5]`. CNOTs from 4 and
5 into 1 compute the parity, resource consumption rotates it, and reversing
those CNOTs restores the original wire labels. The four check wires are finally
measured in X; only `0000` is accepted.

To determine the output phase, let `a` be the four computational bits on the
check wires. Summing `(-1)^(a⋅j)` over the fifteen nonzero binary columns,
with `a⋅j` the binary dot product, gives 15 when `a=0000` and −1 otherwise.
In the product of fifteen rotations,
these phases coincide modulo `2π` on the check register, while the relative
phase of output `|1⟩` to `|0⟩` is `exp(-iπ/4)`. Thus the output is initially
`T†|+⟩`; a final S converts it to `|A⟩`. The circuit simulates that final gate
and its noise too.

```@example factory_reference
reference = distill15()
@assert abs2(dot(target, reference.state)) ≈ 1
@assert reference.raw_states == 15
(attempts=reference.attempts, raw_states=reference.raw_states,
 check_records=reference.syndromes)
```

The six-wire `distill15` function is an unencoded reference circuit. A fault
in its final S can corrupt the output after the checks have passed. Protecting
these Clifford operations requires encoding: `steane_distill15` uses Steane
blocks throughout the distillation circuit and produces a seven-qubit output.

## Acceptance probability and output error

Consider a model in which each of the fifteen raw inputs independently has
a Z defect with probability `p`, and every other operation is perfect. Let
`eⱼ` indicate which raw inputs were defective. The four check outcomes are the
binary column XOR `H e`. An accepted odd-weight pattern adds logical Z to the
output; an accepted even-weight pattern is harmless.

A single bad input cannot pass because every column is nonzero. Two cannot
pass because distinct columns cannot XOR to zero. Three can: the worked
example is columns 1, 2, and 3. Every pair of distinct nonzero columns determines
a unique third column that closes such a triple. There are `15×14/2 = 105`
pairs, and each triple is counted three times, giving 35 undetected triples.
These patterns contribute the leading term `35p³` to the output error.

Writing `q=1-2p`, the exact acceptance and conditional output-error probabilities
for this model are

```math
P_{\mathrm{accept}}=\frac{1+15q^8}{16},\qquad
p_{\mathrm{out}}=\frac{1+15q^8-q^{15}-15q^7}{2(1+15q^8)}
                =35p^3+O(p^4).
```

These are the ideal-Clifford reference formulas of
[Bravyi and Kitaev, Eqs. 35–38](https://arxiv.org/abs/quant-ph/0403025).
The implementation evaluates positive weight-enumerator terms, avoiding
subtraction of nearly equal floating-point numbers for small `p`.

```@example factory_reference
statistics = distillation15_statistics(0.01)
@assert isapprox(statistics.acceptance, 0.8600903336704239)
@assert isapprox(statistics.output_error, 3.60876839653233e-5)
statistics
```

At 1% raw phase error, about 86.0% of attempts pass. The accepted error is about
`3.61×10⁻⁵`, and the mean raw-state cost is `15/P_accept ≈ 17.44`. Rejected
attempts still consume fifteen inputs; this mean cost assumes independent
retries without a cap. These values apply to perfect Clifford gates, readout,
encoded preparation, and recovery. The full factory simulator includes noise
in those operations as well as errors in the raw states.

## Worked example: one and three raw-state errors

```sh
julia --project=. examples/qec/06_distillation_walkthrough.jl
```

This example begins with a noisy physical raw preparation and a
clean run of the six-wire reference. It then inserts one raw Z defect during
the first factory attempt. Its nonzero column label appears in the check
syndrome, so the attempt is discarded. A second, clean attempt passes: the
printed cost is two attempts and 30 raw states.

The next case places Z defects on columns 1, 2, and 3. The check labels cancel,
so the attempt passes, but its odd number of defects changes the output to
`Z|A⟩`. Since `⟨A|Z|A⟩=0`, the output fidelity with the desired state is zero
even though all four checks passed.

The final table evaluates exact ideal-Clifford probabilities for several raw
error rates and compares them with the leading `35p³` approximation. The
physical preparation noise settings and table probabilities are near the top
of the script. The raw-state demonstration includes physical gate and
preparation noise; the table uses the independent raw-Z model. Example 03
below extends the simulation to the encoded factory and its noisy Clifford
operations.

## The encoded resource factory

`steane_distill15` executes the same fifteen rotations on **five Steane blocks**.
The four check blocks are placed first and the output block last, next to the
resource that will be appended. This permutation changes only the storage
order, not the protocol defined above.

Each attempted factory run follows these physical steps:

1. Prepare five logical plus states with `steane_prepare`. Independent blocks
   can be simulated separately until they interact; joining their tensors
   gives their exact tensor product.
2. For each of the fifteen columns, use transversal CNOTs to compute its parity.
   `steane_raw_magic` prepares seven physical zero states, performs H and T on
   the encoder's input wire, applies the noisy physical encoder, and runs noisy
   Shor recovery. `steane_t` consumes that newly prepared resource.
3. Undo the remaining parity and read the four check blocks in logical X using
   physical measurements and classical Hamming decoding. Any nonzero outcome
   discards the whole attempt.
4. Apply the output S correction transversally and run a final noisy Shor round.
   The returned seven-qubit state goes directly into the computational T gadget.

Successive column labels follow Gray-code order: adjacent labels differ in
exactly one bit. One logical CNOT toggles the changed check bit between
successive rotations. Fifteen toggles visit all nonzero columns,
and one final toggle returns to zero. This gives sixteen logical parity CNOTs;
the fifteen Z rotations commute, so their changed order preserves the ideal
distillation operation.

```julia
using BlueTangle
physical_noise = NoiseModel("depolarizing", 1e-5)
resource = steane_distill15(
    noise=physical_noise, prep_error=1e-5, measurement_error=1e-5,
    raw_error=0.01, max_attempts=20,
)
code = steane_code()
data = steane_prepare(:plus; noise=physical_noise,
                     prep_error=1e-5, measurement_error=1e-5)
gate_result = steane_t(code, data.state, resource.state;
                      noise=physical_noise, measurement_error=1e-5)
corrected = shor_qec_round(code, gate_result.state;
                          noise=physical_noise, prep_error=1e-5,
                          measurement_error=1e-5)
readout = steane_measure(code, corrected.state; basis=:X,
                        noise=physical_noise, measurement_error=1e-5)
```

This calculation simulates the noisy preparations, gates, and measurements
in the factory and computation. The returned `resource.state` passes directly
to the logical T gadget. An ideal decoder can be used separately to diagnose
a simulated output.

The largest entangled MPS has 42 physical wires: five data blocks and one
consumed raw block. The **physical peak is 47**, since preparation of that raw
block also uses five Shor ancillas while the other 35 wires are stored. The
small independent preparation registers are simulated separately to avoid
carrying noninteracting tensors through each operation. If another algorithm
register is kept alive during the factory, its wires must be added to this
physical count. Idle noise is not automatically added.

The independent small-register circuits use statevectors. The entangled
factory uses an MPS with singular-value `cutoff=1e-12` by default and no fixed
bond-dimension cap. This is a numerical approximation, separate from the
physical noise model. When changing that model or making accuracy claims,
reduce `cutoff` and check convergence of the quantities being reported;
`cutoff=0.0` retains every numerically nonzero singular value and can cost much
more memory and time.

## Single-fault protection

A single fault in a raw encoder can produce several physical errors in that
one resource block. With otherwise fault-free recovery, the syndrome is removed
and the remaining error acts on its logical qubit. Any resulting logical state
can be expanded in the orthogonal basis `|A⟩, Z|A⟩`. With only that resource
faulty, the second component gives a nonzero outer distillation syndrome.
Acceptance therefore removes that component from the output.

A single fault in a verified preparation or recovery gadget leaves an error
equivalent, up to stabilizers, to at most one physical error per block. During
the subsequent transversal CNOTs and one-qubit gates, a fault can spread between
blocks but stays at that same physical position within each block. Decoded
logical readout tolerates that position's measurement error. The final Shor
round removes an incoming correctable error or, if the fault is in that final
round, may leave a correctable residual. This is the finite single-fault
condition used here.

Two or more physical faults may defeat the inner code, and three defective raw
resources can already produce the accepted logical error described above.
Consequently, the full noisy factory can perform worse than the `35p³`
reference curve. Gate error probabilities, preparation/readout errors, and raw
resource errors must be varied separately when studying that question. The
current examples do not claim a threshold or a fault-tolerant advantage at the
chosen sample size.

An independent five-qubit calculation checks
the Gray-code circuit against the full ideal operator, including T and T†,
and checks both injection outcomes for an arbitrary single faulty logical
resource. The physical reference tests reject all 15 single and 105 double
raw-Z error patterns, and identify all 35 accepted triples with a logical
error. The raw-preparation audit injects all 201 nonidentity Pauli faults at
its preparation, H/T, and encoder locations, then runs the fault-free
Shor circuit and checks that the result is an encoded logical state. The full
factory analysis combines these circuit and algebra checks with the
transversality argument above.

## Acceptance, resource cost, and computation

The encoded factory returns `state`, `attempts`, `raw_states`, `syndromes`,
`history`, and `trace`. Counts include discarded work. An attempt may end by
passing the checks, failing those checks, or exhausting a bounded Shor protocol.
`history` records those outcomes. `raw_states` counts every raw preparation
started, even when a later step aborts. If no attempt succeeds within
`max_attempts`, the factory raises `QECProtocolFailure`.

The factory trace contains the usual physical location plus `factory_attempt`,
`component`, `block`, and `resource` tags. `resource` is the column label 1–15;
the operation's wire indices are local to the named component. `index` is
global across the factory. Child-gadget events additionally carry
`local_index`, the index within that child; direct parity and output-phase
events do not have that field. Fault callbacks should use the component tags
when targeting a particular raw preparation or logical injection. Aborted
child calls retain their location records; outcomes that could not be returned
are `nothing`.

Run the complete factory-to-computation example from the repository root:

```sh
julia --project=. examples/qec/03_noisy_phase_estimation.jl
```

Every trial first runs the encoded factory, then prepares the computational
`|+⟩`, consumes the accepted resource, corrects the data with noisy Shor
extraction, and performs decoded logical measurement. Preparing the
computational block after the factory keeps the physical peak at 47 wires;
keeping that block alive during the factory would raise the count to 54.
Preparation and recovery use physical gate and measurement sequences.

For a perfect T acting on `|+⟩`, logical X and Y have expectations `cos(π/4)`
and `sin(π/4)`. The script estimates these quadratures from fresh noisy trials,
then estimates the phase with `atan(sin_phi, cos_phi)`. Its default two shots
per basis demonstrate the workflow and have large statistical uncertainty.
The settings are physical gate/preparation/readout errors `10⁻⁴`,
extra raw phase error `0.01`, at most five factory attempts, and MPS
`cutoff=1e-12`. Edit those settings at the top of the script.

The reported `factory_attempts` and `raw_states` include discarded work even
when a trial ultimately aborts. Probability estimates use completed shots
only, with aborted shots listed separately. Approximate 95% Wilson intervals
describe the binomial sampling uncertainty of those completed-shot zero
probabilities. Physical-model assumptions and MPS truncation contribute
separate sources of error. The printed ideal-Clifford statistics describe
independent raw Z defects with perfect Clifford operations.

Each fresh physical trial needs a fresh factory output. Reusing a cloned
simulated output for many measurements estimates the measurement probabilities
conditioned on that one factory trajectory. Fresh factory runs also sample
the variation between preparations.

```@docs
shor_code
reed_muller_code
code_distance
raw_magic_state
distillation15_supports
distill15
distillation15_statistics
steane_raw_magic
steane_distill15
```
