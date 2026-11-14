# Measure errors without directly measuring the encoded information. This
# example prepares a Steane block using physical gates, extracts its syndrome
# with verified Shor cat states, and reads out a short noisy logical circuit.
# All simulations are local; no hardware jobs are submitted.
# Run from the repository root:
#   julia --project=. examples/qec/05_noisy_syndrome_cycle.jl

using BlueTangle, LinearAlgebra

# Settings for the stochastic trajectory in part 3. Parts 1 and 2 isolate
# deliberately inserted faults by leaving stochastic noise off.
gate_error = 5e-4
prep_error = 5e-4
measurement_error = 5e-4
noise = NoiseModel("depolarizing", gate_error)
code = steane_code()

# 1. Prepare |+_L> from physical product states and Shor measurements. The
# preparation's initial check outcomes can be random: its recovery selects
# the +1 stabilizer sector while preserving the intended logical state.
prepared = steane_prepare(:plus)
println("Preparation syndrome history: ", prepared.history)
println("Preparation circuit locations: ", length(prepared.trace))

# A physical Y anticommutes with both X and Z on qubit 3. The six syndrome
# bits follow code.generator row order; 1 means the measured check was -1.
# The algebraic prediction below is only a check of the measured result.
incoming_fault = Op("Y", 3)
damaged = apply(prepared.state, incoming_fault)
repaired = shor_qec_round(code, damaged)
expected_syndrome = pauli_syndrome(code, incoming_fault)
@assert repaired.syndrome == expected_syndrome
@assert length(repaired.history) == 2
@assert get_N(repaired.state) == 7

# Decode only to diagnose this ideal-gadget experiment. No decoded state is
# fed back into the computation; extraction simulates the ancilla circuits.
recovered_fidelity = abs2(dot(plus_state(1), code.decode(repaired.state)))
@assert isapprox(recovered_fidelity, 1; atol=1e-10)
println("Incoming Y3: predicted syndrome = ", expected_syndrome)
println("Incoming Y3: measured history = ", repaired.history)
println("Recovery gates: ", [(event.op.name, event.op.qubit)
                            for event in repaired.trace if event.stage == :recovery])
println("Recovered logical |+> fidelity: ", recovered_fidelity)

# 2. Flip one cat readout in the first check of the first round. The callback
# acts immediately before this MZ, so X flips its classical outcome. Fresh
# ancillas occupy wires 8:12; the first four carry the measured cat parity.
# A fault is restricted to an operand of the selected physical location.
readout_fault = location ->
    (location.round == 1 && location.check == 1 &&
     location.kind == :measure && location.stage == :cat_readout &&
     location.op.qubit == code.n + 1) ?
    Op("X", code.n + 1; noisy=false) : Op[]

repeated = shor_qec_round(code, repaired.state; fault=readout_fault)
false_syndrome = [1; zeros(Int, code.m - 1)]
@assert repeated.history == [false_syndrome, zeros(Int, code.m), zeros(Int, code.m)]
@assert isapprox(abs2(dot(plus_state(1), code.decode(repeated.state))), 1; atol=1e-10)
println("One readout fault: measured history = ", repeated.history)
println("Cat preparations per round/check: ", repeated.attempts)

# Two consecutive complete syndromes must agree before recovery. Here the
# first record is false; the second and third agree on no data error. Applying
# recovery after the first record alone would introduce an avoidable fault.
# Each trace entry contains the physical op, kind, stage, round, check,
# attempt, and measurement outcome (nothing at nonmeasurement locations).
first_readout = first(event for event in repeated.trace
                      if event.kind == :measure && event.stage == :cat_readout)
println("Example measured circuit location: ",
        (index=first_readout.index, stage=first_readout.stage,
         round=first_readout.round, check=first_readout.check,
         attempt=first_readout.attempt, gate=first_readout.op.name,
         wire=first_readout.op.qubit, outcome=first_readout.outcome))
println("Extraction locations with the readout fault: ", length(repeated.trace))

# 3. One fully noisy cycle: |0_L> preparation -> transversal logical H ->
# Shor recovery -> destructive logical X readout. Its ideal answer is bit 0,
# meaning eigenvalue +1. The seven physical readout bits are Hamming-decoded
# before their corrected parity is reported as the logical outcome.
#
# Gate depolarization acts independently on each operand of a noisy gate.
# prep_error applies X after a fresh |0>; measurement_error flips an MZ
# outcome. These settings do not add idle faults, leakage, or correlated
# faults at separate locations. No perfect correction is inserted here.
println("Noisy-cycle probabilities (gate/preparation/measurement): ",
        (gate_error, prep_error, measurement_error))
try
    noisy_prepared = steane_prepare(:zero; noise, prep_error, measurement_error)
    after_H = code.apply(noisy_prepared.state, "H", 1; noise)
    corrected = shor_qec_round(code, after_H; noise, prep_error, measurement_error)
    readout = steane_measure(code, corrected.state; basis=:X, noise, measurement_error)
    println("Noisy preparation history: ", noisy_prepared.history)
    println("Noisy recovery history: ", corrected.history)
    println("Noisy recovery cat attempts: ", corrected.attempts)
    println("Physical X-readout bits: ", readout.bits)
    println("Corrected physical bits: ", readout.corrected_bits)
    println("Logical X outcome: ", readout.outcome, " (ideal: 0)")
    println("Recorded preparation/recovery/readout locations: ",
            (length(noisy_prepared.trace), length(corrected.trace), length(readout.trace)))
catch err
    # A heralded protocol failure is an aborted trajectory, not an outcome 0.
    # Report protocol exhaustion here; propagate other exceptions to the caller.
    err isa QECProtocolFailure || rethrow()
    println("Noisy trajectory aborted: ", sprint(showerror, err))
end

# Rerun to sample another trajectory; no RNG seed is imposed. A single shot
# demonstrates the workflow, not a logical error rate or a threshold. With
# one internal fault the gadget may leave a correctable physical error;
# arbitrary multiple faults can corrupt the logical result or cause an abort.
# Example 03 extends these ingredients with a noisy magic-state factory and T.
