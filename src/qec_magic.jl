"""
    steane_readout(bits::AbstractVector)

Decode seven physical Z- or X-measurement bits from one Steane block. Return
`(outcome, bits, corrected_bits, syndrome, corrected_qubit)`. The three Hamming
parity checks locate one erroneous bit; the logical outcome is the parity of
the corrected word (0 for the +1 eigenvalue, 1 for -1). `corrected_qubit=0`
means the word already passed the parity checks.

One physical Pauli before readout or one wrong measurement bit is tolerated.
Two errors can be misdecoded and change the reported logical outcome.
The qubit/check convention is that of `steane_code`.
"""
function steane_readout(bits::AbstractVector)
    length(bits) == 7 || throw(DimensionMismatch("Steane readout requires seven bits"))
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("readout bits must be zero or one"))
    word = Int.(bits)
    checks = stabilizers_to_generator(_STEANE_CHECKS)[1:3, 1:7]
    syndrome = mod.(checks * word, 2)
    qubit = all(iszero, syndrome) ? 0 : findfirst(q -> checks[:, q] == syndrome, 1:7)
    corrected = copy(word)
    qubit != 0 && (corrected[qubit] ⊻= 1)
    return (outcome=mod(sum(corrected), 2), bits=word, corrected_bits=corrected,
            syndrome=syndrome, corrected_qubit=qubit)
end

_qec_join_magic(state::AbstractVector, magic::AbstractVector) = kron(state, magic)
function _qec_join_magic(state::it.MPS, magic::it.MPS)
    # Resources may be copies of a state with the same ITensor index identities.
    # Fresh resource site/link indices keep the two blocks independent.
    resource = it.replace_siteinds(it.sim(it.linkinds, magic), N_MPS(7))
    return it.MPS(vcat(collect(deepcopy(state)), collect(resource)))
end

function _qec_check_steane_register(code, state, block)
    steane_decoder(code) # validate the code's wire and check convention
    get_N(state) == code.n || throw(DimensionMismatch("expected $(code.n) encoded data qubits"))
    1 <= block <= code.k || throw(ArgumentError("block must be between 1 and $(code.k)"))
    magnitude = la.norm(state)
    isfinite(magnitude) && magnitude > 0 || throw(ArgumentError("the input state must have finite positive norm"))
    return nothing
end

# Destructive readout uses physical single-wire measurements. Resetting an
# already measured/disconnected wire below is only bookkeeping before removing
# it from the simulated register; it is not a perfect reset inside the gadget.
function _qec_readout_block(state, wires, trace, context;
                            noise=false, measurement_error=0.0, fault=nothing, cutoff=nothing)
    bits = Int[]
    for q in wires
        state, bit = _qec_run_location(state, Op("MZ", q), context, trace;
                                      kind=:measure, noise=noise,
                                      measurement_error=measurement_error, fault=fault, cutoff=cutoff)
        push!(bits, bit)
        bit == 1 && (state = apply(state, Op("X", q; noisy=false)))
    end
    return state, steane_readout(bits)
end

"""
    steane_measure(code, state; block=1, basis=:Z, noise=false,
                   measurement_error=0.0, fault=nothing, cutoff=nothing)

Destructively measure logical Z or X on one Steane block using seven physical
measurements and classical Hamming decoding. Return a named tuple containing
`state` (all *other* blocks, in their original order), `outcome`, `bits`,
`corrected_bits`, `syndrome`, `corrected_qubit`, and `trace`. Measuring the last
block leaves a zero-qubit statevector of length one.

For X readout, transversal H is applied before Z measurements. `noise` acts on
these basis-change gates. `measurement_error` is an independent bit-flip
probability at each physical measurement. `fault(location)` can return a Pauli
`Op` or a vector of them, restricted to that location's operands; measurement
faults are applied immediately before measurement and gate faults after gates.
`trace` has the same location convention as `shor_qec_round`.
For MPS input, `cutoff` sets the singular-value truncation used in gate
operations; `nothing` keeps their default settings.
It controls numerical accuracy, not physical noise.

The readout tolerates one physical input error OR one circuit/readout fault;
classical processing is assumed reliable. It is not an error-correction step
for the remaining data, and there is no inverse-encoding circuit in this path.
"""
function steane_measure(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS};
                       block::Int=1, basis::Symbol=:Z, noise::Union{Bool,NoiseModel}=false,
                       measurement_error::Real=0.0, fault=nothing, cutoff::Union{Nothing,Real}=nothing)
    _qec_check_steane_register(code, state, block)
    basis in (:X, :Z) || throw(ArgumentError("Steane destructive readout supports basis=:X or :Z"))
    0 <= measurement_error <= 1 || throw(ArgumentError("measurement_error must lie between zero and one"))
    psi = la.normalize(deepcopy(state))
    wires = 7(block-1)+1:7block
    trace = NamedTuple[]
    if basis == :X
        context = (stage=:readout_basis, round=0, check=block, attempt=0)
        for q in wires
            psi, _ = _qec_run_location(psi, Op("H", q), context, trace;
                                       noise=noise, fault=fault, cutoff=cutoff)
        end
    end
    context = (stage=:readout, round=0, check=block, attempt=0)
    psi, result = _qec_readout_block(psi, wires, trace, context;
                                    noise=noise, measurement_error=measurement_error,
                                    fault=fault, cutoff=cutoff)
    remaining = setdiff(collect(1:code.n), wires)
    psi = reduce_and_relabel_qubits(psi, remaining)
    return merge((state=psi,), result, (trace=trace,))
end

"""
    steane_t(code, state, magic; block=1, dagger=false, noise=false,
             measurement_error=0.0, fault=nothing, cutoff=nothing)

Apply logical T by consuming an explicitly supplied seven-qubit encoded magic
state `magic = |A_L>`, where `|A> = (|0> + exp(iπ/4)|1>)/√2`. With
`dagger=true`, supply `|A†_L>` instead and the gadget applies T†.

The circuit appends the resource, applies seven transversal CNOTs with data as
control, destructively measures the resource in Z with Hamming decoding, and
applies logical S (S† for T†) if the decoded outcome is one. It returns
`(state, outcome, bits, corrected_bits, syndrome, corrected_qubit, trace)`;
`state` contains the same data wires as the input, without the consumed resource.
Both arguments must be statevectors or both MPS; MPS evolution stays in tensor
form. Noise and deterministic fault callbacks have the `shor_qec_round` convention.
For MPS input, optional `cutoff` controls singular-value truncation in gate
operations; `nothing` keeps their default settings.
Check convergence when using truncation for quantitative calculations.

This fault-tolerant consumption circuit assumes the supplied magic block has
the intended logical state, up to its specified resource error model. Use
`steane_distill15` to prepare and distill a noisy encoded resource;
`steane_code().encode(T|+>)` supplies an ideal input for simulation checks.
A logical resource error can become a logical data error and is not removed
by stabilizer correction.

With otherwise correct inputs, a single physical fault in a CNOT, readout,
conditional correction, or either input block produces at most a correctable
physical error on the output. Use a fault-tolerant correction round between
gadgets to handle accumulated errors; consumption itself performs no recovery.
"""
function steane_t(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS},
                  magic::Union{AbstractVectorS,it.MPS}; block::Int=1, dagger::Bool=false,
                  noise::Union{Bool,NoiseModel}=false, measurement_error::Real=0.0,
                  fault=nothing, cutoff::Union{Nothing,Real}=nothing)
    _qec_check_steane_register(code, state, block)
    get_N(magic) == 7 || throw(DimensionMismatch("the magic resource must be one seven-qubit Steane block"))
    (state isa it.MPS) == (magic isa it.MPS) || throw(ArgumentError("data and magic must use the same state representation"))
    magnitude = la.norm(magic)
    isfinite(magnitude) && magnitude > 0 || throw(ArgumentError("the magic resource must have finite positive norm"))
    0 <= measurement_error <= 1 || throw(ArgumentError("measurement_error must lie between zero and one"))
    psi = _qec_join_magic(la.normalize(state), la.normalize(magic))
    trace = NamedTuple[]
    context = (stage=:t_cnot, round=0, check=block, attempt=0)
    for q in 1:7
        psi, _ = _qec_run_location(psi, Op("CX", 7(block-1)+q, code.n+q), context, trace;
                                   noise=noise, fault=fault, cutoff=cutoff)
    end
    context = (stage=:t_readout, round=0, check=block, attempt=0)
    psi, result = _qec_readout_block(psi, code.n+1:code.n+7, trace, context;
                                    noise=noise, measurement_error=measurement_error,
                                    fault=fault, cutoff=cutoff)
    psi = reduce_and_relabel_qubits(psi, collect(1:code.n))
    if result.outcome == 1
        context = (stage=:t_correction, round=0, check=block, attempt=0)
        # In this Steane convention physical SD implements logical S.
        for q in 7(block-1)+1:7block
            psi, _ = _qec_run_location(psi, Op(dagger ? "S" : "SD", q), context, trace;
                                       noise=noise, fault=fault, cutoff=cutoff)
        end
    end
    return merge((state=psi,), result, (trace=trace,))
end
