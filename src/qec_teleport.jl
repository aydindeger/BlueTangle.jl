"""
    steane_teleport(code, state; block=1, noise=false, prep_error=0.0,
                    measurement_error=0.0, fault=nothing, max_rounds=6,
                    max_attempts=20, cutoff=1e-12)

Teleport logical qubit `block` of a `steane_code` register onto a fresh
seven-qubit receiver. The receiver replaces the source in the returned wire
order, preserving entanglement with the other blocks. Return
`(state, x_readout, z_readout, trace)`. Each readout is a `steane_readout`
record with the seven physical bits, their Hamming syndrome, the corrected
word and the decoded logical outcome.

`steane_prepare` prepares a sender in logical |+⟩ and a receiver in logical
|0⟩. Transversal CNOT from sender to receiver makes the encoded Bell pair.
A second transversal CNOT from source to sender and destructive source X and
sender Z measurements complete the Bell measurement. The receiver gets
logical X if `z_readout.outcome == 1` and logical Z if `x_readout.outcome == 1`.

With reliable classical processing, a fault-free circuit removes a single
incoming Pauli on the source. With clean encoded input, a Pauli fault at one
physical location can leave at most one erroneous receiver qubit, up to
stabilizers. The Steane code can correct this residual error. Multiple faults
can cause a logical error.

The noise, fault callback and preparation retry limits follow `steane_prepare`.
Preparation exhaustion raises `QECProtocolFailure`. Trace indices count from
one across the entire protocol. Components `:sender_prepare` and
`:receiver_prepare` use local preparation wires. Component `:teleport` uses
the original data wires followed by seven sender and seven receiver wires,
including for feed-forward gates.

Preparation uses independent small statevectors. The joint circuit runs as
an MPS with singular-value truncation `cutoff`. A statevector input returns a
statevector and an MPS input returns an MPS.
"""
function steane_teleport(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS};
                        block::Int=1, noise::Union{Bool,NoiseModel}=false,
                        prep_error::Real=0.0, measurement_error::Real=0.0,
                        fault=nothing, max_rounds::Int=6, max_attempts::Int=20,
                        cutoff::Real=1e-12)
    _qec_check_steane_register(code, state, block)
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative"))
    trace = NamedTuple[]
    psi = la.normalize(state)
    if state isa AbstractVector
        psi = to_MPS(psi, N_MPS(code.n); cutoff)
    end
    for (component, basis) in ((:sender_prepare, :plus), (:receiver_prepare, :zero))
        prepared = _qec_call_with_trace(trace, (; component), fault) do callback
            steane_prepare(basis; noise, prep_error, measurement_error, max_rounds,
                           max_attempts, fault=callback)
        end
        resource = to_MPS(prepared.state, N_MPS(7); cutoff=0.0)
        psi = _qec_join_states(psi, resource)
    end

    source = 7(block-1)+1:7block
    sender = code.n+1:code.n+7
    receiver = code.n+8:code.n+14
    context(stage) = (component=:teleport, stage=stage, round=0, check=block, attempt=0)
    for (stage, controls, targets) in ((:bell_cnot, sender, receiver),
                                       (:teleport_cnot, source, sender))
        for (c, t) in zip(controls, targets)
            psi, _ = _qec_run_location(psi, Op("CX", c, t), context(stage), trace; noise, fault, cutoff)
        end
    end
    for q in source
        psi, _ = _qec_run_location(psi, Op("H", q), context(:readout_basis), trace; noise, fault, cutoff)
    end
    psi, x_readout = _qec_readout_block(psi, source, trace, context(:x_readout);
                                        noise, measurement_error, fault, cutoff)
    psi, z_readout = _qec_readout_block(psi, sender, trace, context(:z_readout);
                                        noise, measurement_error, fault, cutoff)

    for (readout, pauli, stage) in ((z_readout, "X", :x_correction), (x_readout, "Z", :z_correction))
        if readout.outcome == 1
            for q in receiver
                psi, _ = _qec_run_location(psi, Op(pauli, q), context(stage), trace; noise, fault, cutoff)
            end
        end
    end
    # Discard both measured blocks before moving any MPS sites. The receiver
    # now occupies the last seven sites. Relabel it without physical gates.
    retained = vcat(setdiff(1:code.n, source), receiver)
    psi = reduce_and_relabel_qubits(psi, retained)
    if block != code.k
        psi = la.normalize(it.movesites(psi, (code.n-6:code.n) .=> source; cutoff))
    end
    output = state isa it.MPS ? psi : to_state(psi)
    return (state=output, x_readout=x_readout, z_readout=z_readout, trace=trace)
end
