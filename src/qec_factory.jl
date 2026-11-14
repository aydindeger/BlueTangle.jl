"""
    steane_raw_magic(; noise=false, prep_error=0.0, measurement_error=0.0,
                      raw_error=0.0, dagger=false, fault=nothing,
                      max_rounds=6, max_cat_attempts=20)

Prepare a noisy raw encoded |A⟩ (or |A†⟩) from physical |0⟩ preparations,
H and T/TD on the encoder's input wire, the physical Steane encoder, and noisy
Shor recovery. Every preparation, encoding gate, and recovery location uses
the specified physical error model. `raw_error` additionally applies Z just
after the physical T/TD, independently with that probability.

This is deliberately a *raw* resource: an encoder fault can cause a logical
error that Shor recovery cannot remove. Feed it into `steane_distill15` before
using it for a protected computation. Returns the Shor result, including the
complete preparation trace, and `raw_phase_error` indicating the extra Z draw.
"""
function steane_raw_magic(; noise::Union{Bool,NoiseModel}=false,
                          prep_error::Real=0.0, measurement_error::Real=0.0,
                          raw_error::Real=0.0, dagger::Bool=false, fault=nothing,
                          max_rounds::Int=6, max_cat_attempts::Int=20)
    _qec_check_error_probabilities(prep_error, measurement_error)
    isfinite(raw_error) && 0 <= raw_error <= 1 || throw(ArgumentError("raw_error must lie in [0,1]"))
    code = steane_code()
    state = zero_state(7)
    trace = NamedTuple[]
    context = (stage=:raw_prep, round=0, check=0, attempt=0)
    for q in 1:7
        state, _ = _qec_run_location(state, Op("I", q; noisy=false), context, trace;
                                     kind=:prep, prep_error=prep_error, fault=fault)
    end
    q = only(code.logical_indices)
    for name in ("H", dagger ? "TD" : "T")
        context = (stage=name == "H" ? :raw_h : :raw_t, round=0, check=0, attempt=0)
        state, _ = _qec_run_location(state, Op(name, q), context, trace; noise=noise, fault=fault)
    end
    raw_phase_error = rand() < raw_error
    raw_phase_error && (state = apply(state, Op("Z", q; noisy=false)))
    trace[end] = merge(trace[end], (raw_phase_error=raw_phase_error,))
    context = (stage=:raw_encode, round=0, check=0, attempt=0)
    for op in code.ops_encoding
        state, _ = _qec_run_location(state, op, context, trace; noise=noise, fault=fault)
    end
    recovered = _shor_qec_round(code, state, trace; noise=noise, prep_error=prep_error,
        measurement_error=measurement_error, fault=fault,
        max_rounds=max_rounds, max_attempts=max_cat_attempts)
    return merge(recovered, (raw_phase_error=raw_phase_error,))
end

# Expose the exact parity-reuse schedule to independent algebra tests. The
# first fifteen entries precede their column's T injection; column zero is
# the final uncompute. Controls are the four check-block labels, output is 5.
function _steane_distillation15_steps()
    supports = distillation15_supports()
    steps = NamedTuple[]
    previous = Int[]
    for counter in 1:15
        column = xor(counter, counter >> 1)
        current = [q-1 for q in supports[column] if q != 1]
        push!(steps, (column=column, controls=symdiff(previous, current)))
        previous = current
    end
    push!(steps, (column=0, controls=copy(previous)))
    return steps
end

# Child gadgets run on independent small registers where possible. Preserve
# their local wire convention and label the component, block, and raw column.
# Recording before each callback also retains the work of aborted attempts.
function _factory_call(f, trace, tags, fault)
    first = length(trace)
    function record_fault(location)
        event = merge(location, tags,
                      (index=first+location.index, local_index=location.index, outcome=nothing))
        push!(trace, event)
        return fault === nothing ? Op[] : fault(event)
    end
    result = f(record_fault)
    for (i, location) in enumerate(result.trace)
        trace[first+i] = merge(location, tags, (index=first+i, local_index=i))
    end
    return result
end

"""
    steane_distill15(; noise=false, prep_error=0.0, measurement_error=0.0,
                     raw_error=0.0, dagger=false, max_attempts=20,
                     max_rounds=6, max_cat_attempts=20, cutoff=1e-12, fault=nothing)

Run the [[15,1,3]] 15-to-1 distillation circuit on five Steane-encoded blocks,
consuming fifteen explicitly prepared noisy raw encoded resources per complete
attempt. All Clifford gates, preparations, readouts, and final Shor recovery
are physical noisy circuits; no supplied ideal resource or ideal recovery is
used. Returns `(state, attempts, raw_states, syndromes, history, trace)` with a
seven-qubit statevector ready for `steane_t` (`dagger=true` produces |A†⟩).

The four check blocks occupy positions 1:4 and the output is block 5. Nonzero
four-bit column labels are traversed in Gray-code order, reusing the computed
parity: sixteen logical CNOTs implement all fifteen Pauli rotations. Logical
X readout postselects the four checks, and final S (SD for dagger) fixes the
output phase. A failed check or `QECProtocolFailure` restarts the entire attempt;
exhaustion raises `QECProtocolFailure` with the number of raw resources used.

The entangled register is an MPS of at most 42 wires. Independent preparation
circuits use small exact statevectors before tensoring into that register.
The physical footprint is at most 47 wires: 35 stored data wires plus a raw
seven-qubit resource and its five Shor ancillas. Idle noise is not implicit.
`cutoff` controls MPS singular-value truncation during the entangled circuit;
there is no fixed bond-dimension cap. Check convergence when changing the
noise model or circuit. The independent small-register calculations are exact
up to floating-point arithmetic.

`history` includes accepted, rejected, and aborted attempts. `raw_states` counts
all started raw preparations, including discarded/aborted work. The unified
`trace` adds `factory_attempt`, `component`, `block`, and `resource` (column
label 1:15) to each event. Its `op` indices are local to the named component.
Child-gadget events have `local_index` for their child's index; direct parity
and output-phase events do not. `index` counts the whole factory. Use these
labels in `fault` callbacks. Aborted events retain `outcome=nothing`
where a child could not return its completed record.

The encoded circuit contains a single physical fault outside raw preparation
to one position per block. A single faulty raw resource is filtered on
acceptance by the outer distance-three code. This is a single-fault statement,
not a threshold or a guarantee for multiple noisy locations. The familiar
35p³ law applies only to independent raw Z defects and ideal Clifford gadgets;
use `distillation15_statistics` for that reference, not for general gate noise.
"""
function steane_distill15(; noise::Union{Bool,NoiseModel}=false,
                          prep_error::Real=0.0, measurement_error::Real=0.0,
                          raw_error::Real=0.0, dagger::Bool=false,
                          max_attempts::Int=20, max_rounds::Int=6,
                          max_cat_attempts::Int=20, cutoff::Real=1e-12, fault=nothing)
    _qec_check_error_probabilities(prep_error, measurement_error)
    isfinite(raw_error) && 0 <= raw_error <= 1 || throw(ArgumentError("raw_error must lie in [0,1]"))
    max_attempts >= 1 || throw(ArgumentError("max_attempts must be positive"))
    max_rounds >= 2 || throw(ArgumentError("max_rounds must be at least two"))
    max_cat_attempts >= 1 || throw(ArgumentError("max_cat_attempts must be positive"))
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative"))
    code = steane_code(5)
    steps = _steane_distillation15_steps()
    trace, history = NamedTuple[], NamedTuple[]
    syndromes = Vector{Int}[]
    raw_states = 0
    for attempt in 1:max_attempts
        started = raw_states
        try
            # Independent blocks need not share a 35-wire simulation during
            # preparation. Joining their MPS is exactly their tensor product.
            blocks = it.MPS[]
            for b in 1:5
                tags = (factory_attempt=attempt, component=:data_prepare, block=b, resource=0)
                prepared = _factory_call(trace, tags, fault) do callback
                    steane_prepare(:plus; noise=noise, prep_error=prep_error,
                        measurement_error=measurement_error, fault=callback,
                        max_rounds=max_rounds, max_attempts=max_cat_attempts)
                end
                push!(blocks, to_MPS(prepared.state, N_MPS(7); cutoff=0.0))
            end
            state = it.MPS(vcat([collect(block) for block in blocks]...))
            for step in steps[1:end-1]
                column = step.column
                context = (stage=:factory_parity, round=0, check=column, attempt=attempt,
                           factory_attempt=attempt, component=:parity, block=5, resource=column)
                for b in step.controls, op in code.ops("CX", b, 5)
                    state, _ = _qec_run_location(state, op, context, trace; noise=noise, fault=fault, cutoff=cutoff)
                end
                raw_states += 1
                tags = (factory_attempt=attempt, component=:raw, block=0, resource=column)
                raw = _factory_call(trace, tags, fault) do callback
                    steane_raw_magic(; noise=noise, prep_error=prep_error,
                        measurement_error=measurement_error, raw_error=raw_error,
                        dagger=dagger, fault=callback, max_rounds=max_rounds,
                        max_cat_attempts=max_cat_attempts)
                end
                resource = to_MPS(raw.state, N_MPS(7); cutoff=0.0)
                tags = (factory_attempt=attempt, component=:injection, block=5, resource=column)
                injected = _factory_call(trace, tags, fault) do callback
                    steane_t(code, state, resource; block=5, dagger=dagger,
                             noise=noise, measurement_error=measurement_error, fault=callback, cutoff=cutoff)
                end
                state = injected.state
            end
            context = (stage=:factory_parity, round=0, check=0, attempt=attempt,
                       factory_attempt=attempt, component=:parity, block=5, resource=0)
            for b in last(steps).controls, op in code.ops("CX", b, 5)
                state, _ = _qec_run_location(state, op, context, trace; noise=noise, fault=fault, cutoff=cutoff)
            end
            bits = Int[]
            for b in 1:4
                tags = (factory_attempt=attempt, component=:check_readout, block=b, resource=0)
                measured = _factory_call(trace, tags, fault) do callback
                    steane_measure(steane_code(6-b), state; block=1, basis=:X,
                                   noise=noise, measurement_error=measurement_error, fault=callback, cutoff=cutoff)
                end
                state = measured.state
                push!(bits, measured.outcome)
            end
            push!(syndromes, bits)
            if any(!iszero, bits)
                push!(history, (attempt=attempt, accepted=false, reason=:checks,
                                raw_states=raw_states-started, syndrome=bits))
                continue
            end
            # Only seven wires remain. Conversion is exact and avoids using
            # the larger MPS representation for the final physical Shor round.
            output = to_state(state)
            context = (stage=:factory_output_phase, round=0, check=0, attempt=attempt,
                       factory_attempt=attempt, component=:output_phase, block=1, resource=0)
            for op in steane_code().ops(dagger ? "SD" : "S", 1)
                output, _ = _qec_run_location(output, op, context, trace; noise=noise, fault=fault)
            end
            tags = (factory_attempt=attempt, component=:output_correction, block=1, resource=0)
            corrected = _factory_call(trace, tags, fault) do callback
                shor_qec_round(steane_code(), output; noise=noise, prep_error=prep_error,
                    measurement_error=measurement_error, fault=callback,
                    max_rounds=max_rounds, max_attempts=max_cat_attempts)
            end
            push!(history, (attempt=attempt, accepted=true, reason=:accepted,
                            raw_states=raw_states-started, syndrome=bits))
            return (state=corrected.state, attempts=attempt, raw_states=raw_states,
                    syndromes=syndromes, history=history, trace=trace)
        catch err
            err isa QECProtocolFailure || rethrow()
            push!(history, (attempt=attempt, accepted=false, reason=:protocol,
                            raw_states=raw_states-started, syndrome=nothing))
        end
    end
    throw(QECProtocolFailure("Encoded 15-to-1 factory failed after $max_attempts attempts and $raw_states raw preparations"))
end
