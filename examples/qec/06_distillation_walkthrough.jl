# A small, physical 15-to-1 distillation walkthrough.
# Fifteen raw T|+> states are consumed one at a time, so this reference circuit
# needs only six wires. Four X checks reject every single raw Z defect. Three
# suitable defects can pass all checks and leave the output wrong.
# Its Clifford operations are unprotected: use the encoded steane_distill15
# factory in 03_noisy_phase_estimation.jl for the fault-tolerant algorithm.
# Run from the repository root; this script performs local simulation only:
#   julia --project=. examples/qec/06_distillation_walkthrough.jl

using BlueTangle, LinearAlgebra, Printf

# Settings for the first, noisy raw-preparation demonstration. The subsequent
# controlled examples turn stochastic noise off to isolate the filter's action.
gate_error = 0.001
prep_error = 0.001
raw_phase_error = 0.02
reference_probabilities = [1e-4, 1e-3, 1e-2, 0.05]
noise = NoiseModel("depolarizing", gate_error)
A = ComplexF64[1, exp(im*pi/4)] / sqrt(2)       # desired magic state T|+>
ZA = ComplexF64[1, -exp(im*pi/4)] / sqrt(2)     # orthogonal phase-defective state

# 1. Prepare a raw resource using physical |0>, H, and T locations.
# Noise is sampled as a pure-state trajectory. One overlap is therefore one
# realization, not the ensemble fidelity or an estimate of the error rate.
raw = raw_magic_state(; noise, prep_error, phase_error=raw_phase_error)
@assert isapprox(norm(raw.state), 1; atol=1e-12)
@assert [loc.stage for loc in raw.trace] == [:raw_prep, :raw_h, :raw_t]
println("One noisy raw preparation: overlap with A = ", abs2(dot(A, raw.state)))
println("  Extra raw Z sampled: ", raw.trace[end].phase_error_applied)

# 2. With no faults, the injection/readout circuit returns A.
# The simulation prepares and consumes fifteen separate raw states through
# physical gates and measurements.
ideal = distill15(; max_attempts=1)
@assert ideal.attempts == 1 && ideal.raw_states == 15
@assert ideal.syndromes == [zeros(Int, 4)]
@assert isapprox(abs2(dot(A, ideal.state)), 1; atol=1e-12)
println("Ideal circuit: ", ideal.raw_states, " raw states, checks ", ideal.syndromes)

# 3. A single raw Z defect gives a nonzero check string and costs an attempt.
# The callback injects Z after the raw T at column 1 of attempt 1 only. Faults
# target physical locations; all other operations here are noiseless.
one_bad_raw = loc -> loc.stage == :raw_t && loc.attempt == 1 && loc.check == 1 ?
    [Op("Z", loc.op.qubit; noisy=false)] : Op[]
retried = distill15(; fault=one_bad_raw, max_attempts=2)
@assert retried.attempts == 2 && retried.raw_states == 30
@assert retried.syndromes == [[0, 0, 0, 1], zeros(Int, 4)]
@assert count(loc -> loc.stage == :raw_t, retried.trace) == 30
@assert isapprox(abs2(dot(A, retried.state)), 1; atol=1e-12)
println("One raw Z: checks ", retried.syndromes, "; attempts = ", retried.attempts,
        ", raw states consumed = ", retried.raw_states)

# 4. Three raw Z defects can produce an accepted logical error. Check columns
# are the four-bit representations of 1:15. Since 0001 XOR 0010 XOR 0011 = 0000,
# defects on columns 1,2,3 are undetected. Their odd parity leaves a logical Z on A.
three_bad_raw = loc -> loc.stage == :raw_t && loc.check in (1, 2, 3) ?
    [Op("Z", loc.op.qubit; noisy=false)] : Op[]
undetected = distill15(; fault=three_bad_raw, max_attempts=1)
@assert undetected.syndromes == [zeros(Int, 4)]
@assert undetected.attempts == 1 && undetected.raw_states == 15
@assert isapprox(abs2(dot(ZA, undetected.state)), 1; atol=1e-12)
@assert abs2(dot(A, undetected.state)) < 1e-12
println("Three raw Z defects on columns 1,2,3: accepted, but output is Z*A.")

# 5. Exact ideal-Clifford statistics for independent raw Z probability p.
# There are 35 undetected triples. At small p the conditional output error is
# 35p^3 + O(p^4). Rejection also consumes resources: unlimited independent
# retries require 15/P(accept) raw states on average. These exact references
# exclude noisy Cliffords, preparation/readout faults, and finite retry caps.
println("\nIndependent raw-Z reference, with perfect Clifford operations:")
@printf("%9s %12s %14s %14s %12s\n", "p", "P(accept)", "P(error|acc)", "35*p^3", "mean raw")
for p in reference_probabilities
    reference = distillation15_statistics(p)
    @assert 0 < reference.output_error < p
    @assert reference.mean_raw_states >= 15
    @printf("%9.1e %12.8f %14.6e %14.6e %12.5f\n", p,
            reference.acceptance, reference.output_error, 35p^3, reference.mean_raw_states)
end
small_p = first(reference_probabilities)
@assert isapprox(distillation15_statistics(small_p).output_error / (35small_p^3),
                 1; rtol=0.01)

println("\nThis six-wire example explains the filter; it does not protect its Cliffords.")
println("Continue with 03_noisy_phase_estimation.jl for noisy encoded preparation,")
println("the Steane 15-to-1 factory, resource consumption, and Shor correction.")
