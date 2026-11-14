"""
    distillation15_supports()

Return the fifteen Z-product supports of the compressed [[15,1,3]]
Reed–Muller distillation circuit. Wire 1 is the output; wires 2:5 are checks.
Column `j` is `(1, binary(j, 4)...)`, in most-significant-bit order, for
`j=1:15`. Prepare five plus states, apply a T rotation about each Z product,
measure wires 2:5 in X, and accept only four zero outcomes. The ideal output
before the final S is T†|+>; S converts it to T|+>.

The fifteen rotations are implemented by actual magic-state consumption in
`distill15`. This schedule is also shared with the encoded factory.
"""
function distillation15_supports()
    return [[1; [i+1 for i in 1:4 if ((j >> (4-i)) & 1) == 1]] for j in 1:15]
end

"""
    raw_magic_state(; noise=false, prep_error=0.0, phase_error=0.0,
                     dagger=false, fault=nothing)

Physically prepare a noisy single-qubit magic state using |0> preparation,
H, and T (T† for `dagger=true`). `noise` acts after both gates; `prep_error`
is an X fault after preparation. `phase_error` adds an independent Z error
after the T gate, allowing the standard distillation input model to be studied
separately from gate noise. Return `(state, trace)`.

This is a raw, unprotected preparation, not a distilled or encoded resource.
Fault callbacks use the same physical-location convention as `shor_qec_round`.
"""
function raw_magic_state(; noise::Union{Bool,NoiseModel}=false, prep_error::Real=0.0,
                         phase_error::Real=0.0, dagger::Bool=false, fault=nothing)
    _qec_check_error_probabilities(prep_error, phase_error)
    trace = NamedTuple[]
    state = _distillation_raw_magic(zero_state(1), 1, trace, 0, 1;
                                    noise, prep_error, phase_error, dagger, fault)
    return (state=la.normalize(state), trace=trace)
end

function _distillation_raw_magic(state, q, trace, rotation, attempt;
                                 noise, prep_error, phase_error, dagger, fault)
    context = (round=0, check=rotation, attempt=attempt)
    state, _ = _qec_run_location(state, Op("I", q; noisy=false),
                                 merge(context, (stage=:raw_prep,)), trace;
                                 kind=:prep, prep_error, fault)
    state, _ = _qec_run_location(state, Op("H", q),
                                 merge(context, (stage=:raw_h,)), trace; noise, fault)
    state, _ = _qec_run_location(state, Op(dagger ? "TD" : "T", q),
                                 merge(context, (stage=:raw_t,)), trace; noise, fault)
    phase_fault = rand() < phase_error
    phase_fault && (state = apply(state, Op("Z", q; noisy=false)))
    trace[end] = merge(trace[end], (phase_error_applied=phase_fault,))
    return state
end

function _distillation15_attempt(trace, attempt; noise, prep_error, measurement_error,
                                 raw_error, dagger, fault)
    state = zero_state(5)
    context = (round=0, check=0, attempt=attempt)
    for q in 1:5
        state, _ = _qec_run_location(state, Op("I", q; noisy=false),
                                     merge(context, (stage=:distill_prep,)), trace;
                                     kind=:prep, prep_error, fault)
        state, _ = _qec_run_location(state, Op("H", q),
                                     merge(context, (stage=:distill_basis,)), trace; noise, fault)
    end
    for (rotation, support) in enumerate(distillation15_supports())
        context = (round=0, check=rotation, attempt=attempt)
        # Compute this Z product onto wire 1. After uncomputing, a T on wire 1
        # becomes exp(-i*pi*Z_product/8), up to a global phase.
        for q in support[2:end]
            state, _ = _qec_run_location(state, Op("CX", q, 1),
                                         merge(context, (stage=:distill_compute,)), trace; noise, fault)
        end
        state = kron(state, zero_state(1))
        state = _distillation_raw_magic(state, 6, trace, rotation, attempt;
                                        noise, prep_error, phase_error=raw_error, dagger, fault)
        state, _ = _qec_run_location(state, Op("CX", 1, 6),
                                     merge(context, (stage=:distill_inject,)), trace; noise, fault)
        state, bit = _qec_run_location(state, Op("MZ", 6),
                                       merge(context, (stage=:distill_inject_readout,)), trace;
                                       kind=:measure, measurement_error, fault)
        if bit == 1
            state, _ = _qec_run_location(state, Op(dagger ? "SD" : "S", 1),
                                         merge(context, (stage=:distill_feedforward,)), trace; noise, fault)
            state = apply(state, Op("X", 6; noisy=false))
        end
        state = reduce_and_relabel_qubits(state, collect(1:5))
        for q in reverse(support[2:end])
            state, _ = _qec_run_location(state, Op("CX", q, 1),
                                         merge(context, (stage=:distill_uncompute,)), trace; noise, fault)
        end
    end
    bits = Int[]
    context = (round=0, check=0, attempt=attempt)
    for q in 2:5
        state, _ = _qec_run_location(state, Op("H", q),
                                     merge(context, (stage=:distill_check_basis,)), trace; noise, fault)
        state, bit = _qec_run_location(state, Op("MZ", q),
                                       merge(context, (stage=:distill_check_readout,)), trace;
                                       kind=:measure, measurement_error, fault)
        push!(bits, bit)
        bit == 1 && (state = apply(state, Op("X", q; noisy=false)))
    end
    if any(!iszero, bits)
        return (state=nothing, accepted=false, syndrome=bits)
    end
    state = reduce_and_relabel_qubits(state, [1])
    # All fifteen T-product rotations produce T† on the output. This last
    # physical Clifford is intentionally noisy too; its faults are not caught
    # by earlier checks in this unencoded reference circuit.
    state, _ = _qec_run_location(state, Op(dagger ? "SD" : "S", 1),
                                 merge(context, (stage=:distill_output,)), trace; noise, fault)
    return (state=la.normalize(state), accepted=true, syndrome=bits)
end

"""
    distill15(; noise=false, prep_error=0.0, measurement_error=0.0,
               raw_error=0.0, dagger=false, max_attempts=20, fault=nothing)

Run the compressed 15-to-1 Reed–Muller magic-state distillation circuit with
physical noisy operations. Each attempt prepares five plus-state data wires
and consumes fifteen newly prepared raw magic states, one at a time on a sixth
wire. Four X measurements decide acceptance. Rejected attempts are discarded
and repeated; exhaustion raises `QECProtocolFailure`.

Return `(state, attempts, raw_states, syndromes, trace)`, where `state` is the
accepted single-qubit T|+> resource (T†|+> for `dagger=true`). Counts and trace
include rejected attempts. Preparation, H/T, parity CNOTs, injection readout,
feed-forward, output S, and final check readout all have explicit noisy
locations. `raw_error` adds independent Z defects after raw T preparations.
Fault callbacks have the Shor convention; `check` identifies rotation 1:15.

This six-wire circuit is an unencoded physical reference. Its Clifford faults
are not protected, so its output is not guaranteed to improve when gate noise
is appreciable. In particular, the ideal-Clifford `35p^3` law does not describe
the full noisy circuit. Use the encoded factory for protected Clifford gates.
"""
function distill15(; noise::Union{Bool,NoiseModel}=false, prep_error::Real=0.0,
                   measurement_error::Real=0.0, raw_error::Real=0.0,
                   dagger::Bool=false, max_attempts::Int=20, fault=nothing)
    _qec_check_error_probabilities(prep_error, measurement_error)
    _qec_check_error_probabilities(raw_error, 0.0)
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive"))
    trace, syndromes = NamedTuple[], Vector{Int}[]
    for attempt in 1:max_attempts
        result = _distillation15_attempt(trace, attempt;
                                         noise, prep_error, measurement_error, raw_error, dagger, fault)
        push!(syndromes, result.syndrome)
        if result.accepted
            return (state=result.state, attempts=attempt, raw_states=15attempt,
                    syndromes=syndromes, trace=trace)
        end
    end
    throw(QECProtocolFailure("15-to-1 distillation rejected $max_attempts attempts"))
end

"""
    distillation15_statistics(p)

Exact reference statistics for independent Z defects of probability `p` on
fifteen raw magic inputs and otherwise perfect Clifford operations. Return
`(acceptance, output_error, mean_raw_states)`. This enumerator is a comparison
for the circuit simulator, never a replacement for its noisy operations.

An input error subset is accepted when its four-bit column XOR vanishes.
An accepted subset of odd size produces a logical Z defect; the smallest such
subsets have size three, and there are 35 of them. Thus output error is
`35p^3 + O(p^4)` near zero, conditioned on acceptance.
"""
function distillation15_statistics(p::Real)
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must lie in [0,1]"))
    # Weight enumerators of the accepted even and odd subsets. The explicit
    # polynomial avoids cancellation of nearly equal terms when p is small.
    correct = sum(count * p^weight * (1-p)^(15-weight)
                  for (weight, count) in ((0,1), (4,105), (6,280), (8,435), (10,168), (12,35)))
    incorrect = sum(count * p^weight * (1-p)^(15-weight)
                    for (weight, count) in ((3,35), (5,168), (7,435), (9,280), (11,105), (15,1)))
    acceptance = correct + incorrect
    return (acceptance=acceptance, output_error=incorrect/acceptance,
            mean_raw_states=15/acceptance)
end
