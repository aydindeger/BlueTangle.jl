# Circuit locations shared by the Shor and magic-state gadgets. Preparation is
# a fresh |0>, with an optional X preparation error. Measurement error is an X
# just before MZ; this is equivalent to a classical bit flip for the data branch.
"""A heralded failure to obtain verified ancillas or a repeated syndrome within the requested bounds."""
struct QECProtocolFailure <: Exception
    message::String
end
Base.showerror(io::IO, failure::QECProtocolFailure) = print(io, failure.message)

function _qec_run_location(state, op::Op, context::NamedTuple, record::Vector;
                           kind::Symbol=:gate, noise::Union{Bool,NoiseModel}=false,
                           prep_error::Real=0.0, measurement_error::Real=0.0,
                           fault=nothing, offset::Int=0, cutoff::Union{Nothing,Real}=nothing)
    kind in (:gate, :prep, :measure) || throw(ArgumentError("unknown QEC location kind: $kind"))
    cutoff === nothing || isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative"))
    tensor_kwargs = state isa it.MPS && cutoff !== nothing ? (; cutoff=cutoff) : (;)
    location = merge(context, (index=length(record)+1, kind=kind, op=op))
    # Cat preparation can run on its five local wires before being joined to
    # the data, while the recorded schedule retains absolute physical indices.
    local_op(o) = offset == 0 ? o : Op(o.name, o.mat, o.qubit-offset,
                                      o.target_qubit < 0 ? -1 : o.target_qubit-offset;
                                      noisy=o.noisy)
    function inject(psi)
        fault === nothing && return psi
        errors = fault(location)
        errors = errors isa Op ? [errors] : errors
        errors isa AbstractVector || throw(ArgumentError("fault callback must return a Pauli Op or vector of them"))
        _qec_pauli_vector(errors, get_N(psi)+offset)
        operands = op.q == 1 ? (op.qubit,) : (op.qubit, op.target_qubit)
        all(e -> e.qubit in operands, errors) ||
            throw(ArgumentError("a location fault may only act on that location's operands"))
        for e in errors
            psi = apply(psi, local_op(e); tensor_kwargs...)
        end
        return psi
    end
    bit = nothing
    if kind == :gate
        state = apply(state, local_op(op); noise=noise, tensor_kwargs...)
        state = inject(state)
    elseif kind == :prep
        rand() < prep_error && (state = apply(state, Op("X", op.qubit-offset; noisy=false)))
        state = inject(state)
    else
        op.name == "MZ" || throw(ArgumentError("QEC measurement locations must use MZ"))
        rand() < measurement_error && (state = apply(state, Op("X", op.qubit-offset; noisy=false)))
        state = inject(state)
        state, measured = apply(state, local_op(op); track_measurements=true)
        bit = only(measured)
    end
    push!(record, merge(location, (outcome=bit,)))
    return state, bit
end

function _qec_check_error_probabilities(prep_error, measurement_error)
    all(p -> isfinite(p) && 0 <= p <= 1, (prep_error, measurement_error)) ||
        throw(ArgumentError("preparation and measurement error probabilities must lie in [0,1]"))
end

"""
    shor_syndrome_circuit(code::StabilizerCode, check::Int)

Return the physical locations for one verified Shor measurement of a Steane
weight-four check. Each entry is `(kind, stage, op)`. The five ancillary wires
are `code.n+1:code.n+5`; `:prep` entries denote freshly prepared |0> states.

The cat uses H(1), CX(1,2), CX(1,3), CX(1,4), in that order. A fifth qubit
verifies Z(1)Z(2) before the cat touches the data. Reject outcome one and restart
the ancilla preparation. This particular verifier is sufficient for this fixed
preparation order; changing that order requires a new fault analysis. Each cat
qubit subsequently interacts with just one data qubit. H and MZ read out the
cat's X parity. The verification measurement is not part of that parity.

This is the schedule used by `shor_qec_round`, including preparation and
measurement locations. It is not a flat circuit that can be run without the
verification rejection and repeated-syndrome control flow.
"""
function shor_syndrome_circuit(code::StabilizerCode, check::Int)
    steane_decoder(code) # Validate the generator and contiguous block layout.
    1 <= check <= code.m || throw(ArgumentError("check must lie between 1 and $(code.m)"))
    row = code.generator[check, :]
    support = [q for q in 1:code.n if row[q] != 0 || row[code.n+q] != 0]
    length(support) == 4 || throw(ArgumentError("Shor extraction requires a weight-four Steane check"))
    a, v = code.n .+ collect(1:4), code.n+5
    circuit = NamedTuple[]
    for q in code.n+1:code.n+5
        push!(circuit, (kind=:prep, stage=:cat_prep, op=Op("I", q; noisy=false)))
    end
    for op in [Op("H", a[1]); [Op("CX", a[1], a[j]) for j in 2:4]]
        push!(circuit, (kind=:gate, stage=:cat_encode, op=op))
    end
    for op in (Op("CX", a[1], v), Op("CX", a[2], v))
        push!(circuit, (kind=:gate, stage=:cat_verify, op=op))
    end
    push!(circuit, (kind=:measure, stage=:cat_verify, op=Op("MZ", v)))
    for (j, q) in enumerate(support)
        name = row[q] == 1 ? "CX" : "CZ"
        push!(circuit, (kind=:gate, stage=:coupling, op=Op(name, a[j], q)))
    end
    for q in a
        push!(circuit, (kind=:gate, stage=:cat_readout, op=Op("H", q)))
    end
    for q in a
        push!(circuit, (kind=:measure, stage=:cat_readout, op=Op("MZ", q)))
    end
    return circuit
end

_qec_join_states(data::AbstractVector, ancilla::AbstractVector) = kron(data, ancilla)
_qec_join_states(data::it.MPS, ancilla::it.MPS) = it.MPS(vcat(collect(deepcopy(data)), collect(deepcopy(ancilla))))
_qec_cat_zero(::AbstractVector) = zero_state(5)
_qec_cat_zero(::it.MPS) = zero_state(N_MPS(5))

function _qec_shor_check(code, state, circuit, round, check, record;
                         noise, prep_error, measurement_error, max_attempts, fault)
    for attempt in 1:max_attempts
        context = (round=round, check=check, attempt=attempt)
        cat = _qec_cat_zero(state)
        verification = 1
        # Only the cat and verifier participate before entry 13. Simulating
        # their preparation separately saves a factor 2^n without changing it.
        for entry in circuit[1:12]
            cat, bit = _qec_run_location(cat, entry.op, merge(context, (stage=entry.stage,)), record;
                                        kind=entry.kind, noise=noise, prep_error=prep_error,
                                        measurement_error=measurement_error, fault=fault, offset=code.n)
            entry.kind == :measure && (verification = bit)
        end
        verification == 1 && continue
        joint = _qec_join_states(state, cat)
        parity = 0
        for entry in circuit[13:end]
            joint, bit = _qec_run_location(joint, entry.op, merge(context, (stage=entry.stage,)), record;
                                          kind=entry.kind, noise=noise, prep_error=prep_error,
                                          measurement_error=measurement_error, fault=fault)
            if entry.kind == :measure
                parity ⊻= bit
                # The measured wire is now a known product factor. Mapping it
                # to zero here is bookkeeping for discarding it, not a gate on
                # data or an unmodelled physical reset used later by the gadget.
                bit == 1 && (joint = apply(joint, Op("X", entry.op.qubit; noisy=false)))
            end
        end
        return reduce_and_relabel_qubits(joint, collect(1:code.n)), parity, attempt
    end
    throw(QECProtocolFailure("Shor cat verification failed after $max_attempts attempts (round $round, check $check)"))
end

"""
    shor_qec_round(code, state; noise=false, prep_error=0.0,
                   measurement_error=0.0, decoder=nothing, max_rounds=6,
                   max_attempts=20, fault=nothing)

Run a noisy, distance-three Shor error-correction gadget for independent Steane
blocks. Fresh verified four-qubit cat states measure the checks; complete
syndromes are repeated until two consecutive records agree. The matching record
is decoded and physical Pauli recovery gates receive the supplied gate noise.
The result is `(state, syndrome, history, attempts, trace)`, with only the
original data wires retained. `attempts[r][j]` counts cat preparations for check
`j` in round `r`. Exhausting either bound raises `QECProtocolFailure`, never a
success flag.

`noise::NoiseModel` applies to actual one- and two-qubit gates, including cat
preparation, verification, data interactions, readout H gates, and recovery.
`prep_error` is the independent probability of X after a fresh |0> preparation;
`measurement_error` is an independent outcome flip at each MZ. No idle faults,
leakage, or correlated faults at distinct locations are introduced implicitly.
Statevectors and native MPS trajectories are supported; density matrices are not.

For reproducible adversarial tests, `fault(location)` returns a Pauli `Op` or
vector of them on that location's operands. A location contains `index`, `kind`,
`op`, `stage`, `round`, `check`, and `attempt`. Faults occur after preparations
and gates, and before measurements. `trace` records these events plus the
measurement `outcome` (or `nothing`). Setting noise and error probabilities to
zero isolates the injected fault. The callback is also invoked on recovery.

The single-fault claim concerns this fixed circuit and repetition rule: with
clean encoded input, one faulty location leaves an error equivalent, up to
stabilizers, to at most one Pauli per block; without internal faults it corrects
one input Pauli per block. It
does not assert success for arbitrary multiple faults or establish a threshold.
Unlike `qec_round`, this function does not replace extraction by ideal check
projectors. The fault-containment claim uses the default `steane_decoder(code)`;
syndrome consistency alone cannot certify the logical action of a custom decoder.
"""
function shor_qec_round(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS}; kwargs...)
    return _shor_qec_round(code, state, NamedTuple[]; kwargs...)
end

function _shor_qec_round(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS}, record::Vector;
                         noise::Union{Bool,NoiseModel}=false, prep_error::Real=0.0,
                         measurement_error::Real=0.0, decoder=nothing,
                         max_rounds::Int=6, max_attempts::Int=20, fault=nothing)
    get_N(state) == code.n || throw(DimensionMismatch("Shor recovery requires $(code.n) data qubits"))
    _qec_check_error_probabilities(prep_error, measurement_error)
    max_rounds >= 2 || throw(ArgumentError("max_rounds must allow at least two syndrome rounds"))
    max_attempts >= 1 || throw(ArgumentError("max_attempts must be positive"))
    state_norm = la.norm(state)
    isfinite(state_norm) && state_norm > 0 || throw(ArgumentError("the input state must have finite positive norm"))
    state = la.normalize(state)
    default_decoder = steane_decoder(code)
    decoder = decoder === nothing ? default_decoder : decoder
    circuits = [shor_syndrome_circuit(code, check) for check in 1:code.m]
    history, attempts = Vector{Int}[], Vector{Int}[]
    for round in 1:max_rounds
        bits, counts = Int[], Int[]
        for check in 1:code.m
            state, bit, count = _qec_shor_check(code, state, circuits[check], round, check, record;
                                              noise=noise, prep_error=prep_error,
                                              measurement_error=measurement_error,
                                              max_attempts=max_attempts, fault=fault)
            push!(bits, bit)
            push!(counts, count)
        end
        push!(history, bits)
        push!(attempts, counts)
        if round >= 2 && bits == history[end-1]
            recovery = _qec_recovery(bits, code.generator, decoder)
            recovery isa AbstractString && throw(ArgumentError(recovery))
            for op in (recovery isa Op ? [recovery] : recovery)
                op.name == "I" && continue
                physical = Op(op.name, op.qubit; noisy=true)
                context = (stage=:recovery, round=round, check=0, attempt=0)
                state, _ = _qec_run_location(state, physical, context, record;
                                              noise=noise, fault=fault)
            end
            return (state=state, syndrome=bits, history=history,
                    attempts=attempts, trace=record)
        end
    end
    throw(QECProtocolFailure("Shor syndrome did not repeat in $max_rounds rounds"))
end
