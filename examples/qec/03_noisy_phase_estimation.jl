# Full noisy cycle: physical raw magic preparation -> encoded 15-to-1 factory
# -> logical T consumption -> Shor correction -> decoded logical readout.
# This Ramsey/phase-estimation experiment measures X and Y after T|+>, whose
# ideal means are cos(pi/4) and sin(pi/4). Every trial makes a fresh resource.
# No ideal encoded magic input, inverse encoder, or ideal QEC is used below.
# Run from the repository root:
#   julia --project=. examples/qec/03_noisy_phase_estimation.jl

using BlueTangle, LinearAlgebra

# Two shots per basis keep the complete factory demonstration manageable.
# Increase this number for meaningful phase statistics; first-run compilation
# and the 42-wire MPS factory make this more expensive than the earlier example.
shots = 2
gate_error = 0.0001                 # per operand, including two-qubit gates
prep_error = 0.0001                 # X after each physical |0> preparation
measurement_error = 0.0001          # flip of each physical measurement outcome
raw_error = 0.01                    # extra Z after each raw physical T
factory_max_attempts = 5
cutoff = 1e-12                      # MPS singular-value truncation; no bond cap
noise = NoiseModel("depolarizing", gate_error)
code = steane_code()

# Verify code properties from actual checks, rather than trusting stored d.
rm_distance = code_distance(reed_muller_code(); max_weight=7)
@assert (rm_distance.distance, rm_distance.dx, rm_distance.dz) == (3, 7, 3)
println("Reed-Muller distance search: d=3, dX=7, dZ=3")

# This reference excludes all physical gate/preparation/readout errors. It is
# useful for understanding the distillation filter, not as a prediction for
# the full noisy factory simulated below.
println("Ideal-Clifford, independent-raw-Z reference: ", distillation15_statistics(raw_error))

function phase_trial(code, basis; noise, prep_error, measurement_error,
                     raw_error, factory_max_attempts, cutoff, count_factory)
    # Build the factory first. Keeping algorithm data live throughout would
    # increase the physical peak from 47 to 54 wires, even without idle noise.
    factory = steane_distill15(; noise, prep_error, measurement_error, raw_error,
                               max_attempts=factory_max_attempts, cutoff,
                               fault=count_factory)
    prepared = steane_prepare(:plus; noise, prep_error, measurement_error)
    injected = steane_t(code, prepared.state, factory.state; noise, measurement_error)
    corrected = shor_qec_round(code, injected.state; noise, prep_error, measurement_error)
    # S† before an X measurement measures the original Y quadrature.
    state = basis == :Y ? code.apply(corrected.state, "SD", 1; noise) : corrected.state
    measured = steane_measure(code, state; basis=:X, noise, measurement_error)
    return (bit=measured.outcome, factory_attempts=factory.attempts,
            raw_states=factory.raw_states, factory_history=factory.history)
end

results = Dict{Symbol,Any}()
for basis in (:X, :Y)
    bits = Int[]
    aborted = total_attempts = total_raw = 0
    for shot in 1:shots
        # These counters include work in a factory that ultimately aborts.
        # The callback returns no injected fault and leaves stochastic noise on.
        attempt_count, raw_count = Ref(0), Ref(0)
        count_factory = loc -> begin
            attempt_count[] = max(attempt_count[], loc.factory_attempt)
            if loc.component == :raw && loc.stage == :raw_prep && loc.op.qubit == 1
                raw_count[] += 1
            end
            Op[]
        end
        println("Running ", basis, " trial ", shot, "/", shots, " with a fresh factory...")
        flush(stdout)
        try
            result = phase_trial(code, basis; noise, prep_error, measurement_error,
                                 raw_error, factory_max_attempts, cutoff, count_factory)
            @assert result.raw_states == raw_count[]
            @assert result.factory_attempts == attempt_count[]
            push!(bits, result.bit)
            println("  logical bit=", result.bit, ", raw resources=", result.raw_states,
                    ", factory attempts=", result.factory_attempts)
            println("  factory history=", [(attempt=r.attempt, reason=r.reason,
                                            check_bits=r.syndrome) for r in result.factory_history])
        catch err
            err isa QECProtocolFailure || rethrow()
            aborted += 1
            println("  heralded failure: ", sprint(showerror, err))
        end
        total_attempts += attempt_count[]
        total_raw += raw_count[]
    end
    completed = length(bits)
    p_zero = completed == 0 ? NaN : count(iszero, bits) / completed
    # A plug-in standard error would be zero if two shots agree. The Wilson
    # score interval remains wide in that case and exposes the limited sample.
    # This approximate 95% binomial interval concerns completed shots only.
    z95 = 1.959963984540054
    interval95 = if completed == 0
        (NaN, NaN)
    else
        denominator = 1 + z95^2/completed
        center = (p_zero + z95^2/(2completed)) / denominator
        radius = z95 * sqrt(p_zero*(1-p_zero)/completed + z95^2/(4completed^2)) / denominator
        (max(0.0, center-radius), min(1.0, center+radius))
    end
    results[basis] = (; p_zero, wilson_interval95=interval95, completed,
                       aborted, factory_attempts=total_attempts, raw_states=total_raw)
    println("Logical ", basis, " results: ", results[basis])
end

cos_phi = 2results[:X].p_zero - 1
sin_phi = 2results[:Y].p_zero - 1
phase_estimate = cos_phi == 0 && sin_phi == 0 ? NaN : atan(sin_phi, cos_phi)
println("Estimated phase (radians): ", phase_estimate, "; ideal: ", pi/4)
println("Ideal probability of zero in either basis: ", (1+inv(sqrt(2)))/2)
println("Probability and phase estimates are conditioned on completed shots.")
println("Wilson intervals describe finite-shot uncertainty, not physical or MPS bias.")
println("Factory/raw counts include rejected attempts and aborted trials.")

# This small run demonstrates the complete physical workflow and its records.
# It does not estimate a threshold or show a noisy advantage. A bad accepted
# resource is possible when multiple faults occur; postselection is not a
# fidelity certificate. Idle noise, leakage, and classical-control faults are
# outside this model. Check cutoff convergence before quantitative MPS studies.
