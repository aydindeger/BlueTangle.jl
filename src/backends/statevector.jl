# Statevector trajectories.

"""
`apply(state::AbstractVectorS, op::QuantumOps)`

Return the state after applying an operation. Measurements and noise sample one
trajectory. Pass `track_measurements=true` to return `(state, bits)`.
"""
function apply(
    state::AbstractVectorS,
    op::QuantumOps;
    noise::Union{NoiseModel,Bool}=false,
    track_measurements::Bool=false,
    kwargs...
)

    _unsupported_keyword_check(kwargs,[:pars]) #:pars are for variational gates

    N=get_N(state)
    mid_measurements=Int[]

    # if op.q!=1 && abs(op.qubit-op.target_qubit)>1
    #     throw("non-local gate $(op.name) is not allowed!")
    # end

    if isa(op,OpF)
        if op.data isa Vector{<:QuantumOps}
            if track_measurements
                state, mid_measurements = apply(state, op.data; noise=noise, track_measurements=true, kwargs...)
            else
                state = apply(state, op.data; noise=noise, kwargs...)
            end
        else
            state = op.apply(state)
        end
    elseif isa(op,OpQC)
        state=op.apply(state)
    elseif op.type=="🔬"
        if isa(op,ifOp)
            state,record=op.born_apply(state,noise;track_measurements,kwargs...)
            track_measurements && append!(mid_measurements,record)
        else
            state,ind=_born_measure(state,op)
            track_measurements && push!(mid_measurements,ind)
        end
        # println("measurement result=$(ind)")
    elseif op.type=="f" #for variational circuits
        state=op.expand(N,kwargs[:pars]...)*state
    else #good old gates
        state=op.expand(N)*state
    end

    ##apply noise.
    if isa(noise, NoiseModel)
        state=apply_noise(state,op,noise)
    end

    return track_measurements ? (state,mid_measurements) : state

end


function _born_measure(state::AbstractVectorS,o::QuantumOps)

    N=get_N(state)
    rotMat=o.expand(N)
    state=rotMat*state#rotate
    state,ind=born_measure_Z(N,state,o.qubit)
    state=rotMat'*state#rotate back

    return state,ind

end


"""
    born_measure_Z(N::Int, state::AbstractVectorS, qubit::Int)

Sample a Z measurement of one qubit in a normalized `N`-qubit statevector.
Return `(collapsed_state, bit)`, with bit 0 for `|0⟩` and 1 for `|1⟩`.
The returned state is normalized and retains all qubits. Qubit 1 is the most
significant tensor factor. For other bases or automatic record collection, use
`apply(state, Op("MZ", qubit); track_measurements=true)` and its MX/MY variants.
"""
function born_measure_Z(N::Int,state::AbstractVectorS,qubit::Int)

    born_ops=[gate.P0, gate.P1]

    #according to benchmark
    # if N<=12
        # prob0=sum(abs2.(hilbert(N,born_ops[1],qubit)*state))
    # else
        prob0=real(BlueTangle.partial_trace(state,qubit))[1,1]
    # end

    ind=rand() < prob0 ? 0 : 1
    return sa.normalize(hilbert(N,born_ops[ind+1],qubit)*state),ind

end


"""
    born_measure_Z(N::Int,state::AbstractVectorS,qubit1::Int,qubit2::Int)

Sample a joint Z measurement on two qubits of a normalized statevector and
return `(collapsed_state, outcome_index)`. This overload's second return value
is a **1-based index**, ordered as 00, 10, 01, 11 for `(qubit1, qubit2)`.
It is not the 0/1 bit returned by the single-qubit overload.
"""
function born_measure_Z(N::Int,state::AbstractVectorS,qubit1::Int,qubit2::Int)

    slist=[expand_multi_op("P0,P0",[qubit1,qubit2],N)*state
    ,expand_multi_op("P1,P0",[qubit1,qubit2],N)*state
    ,expand_multi_op("P0,P1",[qubit1,qubit2],N)*state
    ,expand_multi_op("P1,P1",[qubit1,qubit2],N)*state]

    measure0i0=abs(state'slist[1])
    measure1i0=abs(state'slist[2])
    measure0i1=abs(state'slist[3])
    measure1i1=abs(state'slist[4])

    probs=[measure0i0,measure1i0,measure0i1,measure1i1]
    ind=BlueTangle._weighted_sample(probs)

    return sa.normalize(slist[ind]),ind
end


function _reset_Z(state::AbstractVectorS,qubit::Int)
    N=get_N(state)
    state,ind=born_measure_Z(N,state,qubit)
    if ind==1
        state=Op("X",qubit)*state
    end
    return state,ind
end


##
##========== state preparation ==========

function _state_from_backend(backend, sv_fn, mps_fn, camps_fn)
    backend_name = _resolve_backend(backend)
    if backend_name == "statevector"
        return sv_fn()
    elseif backend_name == "density_matrix"
        s = sv_fn()
        return s * s'
    elseif backend_name == "tensor"
        return mps_fn()
    elseif backend_name == "camps"
        return camps_fn()
    elseif backend_name == "phasepoly"
        throw(ArgumentError("PhasePoly represents circuits directly; use `apply(ops; backend=:phasepoly)` or `phasepoly(ops)`"))
    end
end

"""
    product_state(bits::Vector; backend=:statevector)

Return the computational-basis product state specified by 0/1 `bits` (qubit 1
is the first entry / most significant bit) in the requested `backend`
(`:statevector`, `:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
function product_state(list_of_qubits::Vector; backend::Union{Symbol,AbstractString}=:statevector)
    bits = Int.(sign.(list_of_qubits))
    return _state_from_backend(
        backend,
        () -> sa.sparse(foldl(kron, _bin2state.(bits))),
        () -> product_state(N_MPS(length(bits)), bits),
        () -> begin
            st = CAMPSState(length(bits))
            for (q, b) in enumerate(bits)
                b == 1 && (st = apply(st, Op("X", q)))
            end
            st
        end,
    )
end

"""
    neel_state01(N::Int; backend=:statevector)

Return the `N`-qubit Néel state `|0101…⟩` in the requested `backend`
(`:statevector`, `:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
neel_state01(N::Int; backend::Union{Symbol,AbstractString}=:statevector) =
    product_state([isodd(i) ? 0 : 1 for i=1:N]; backend=backend)

"""
    neel_state10(N::Int; backend=:statevector)

Return the `N`-qubit Néel state `|1010…⟩` in the requested `backend`
(`:statevector`, `:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
neel_state10(N::Int; backend::Union{Symbol,AbstractString}=:statevector) =
    product_state([isodd(i) ? 1 : 0 for i=1:N]; backend=backend)


"""
    Z3(N::Int, j=1; backend=:statevector)
"""
Z3(N::Int, j=1; backend::Union{Symbol,AbstractString}=:statevector) =
    product_state([mod(i+2-j+1,3)==0 ? 1 : 0 for i=1:N]; backend=backend)

"""
    Z3_s(N::Int)
"""
Z3_s(N::Int)=la.normalize(Z3(N,1)+Z3(N,2)+Z3(N,3))

"""
    neel_state_s(N::Int)=(neel_state01(N)+neel_state10(N))/sqrt(2)
"""
neel_state_s(N::Int)=(neel_state01(N)+neel_state10(N))/sqrt(2)

"""
    random_product_state(N::Int; backend=:statevector)

Return a random computational-basis product state on `N` qubits in the
requested `backend` (`:statevector`, `:density_matrix`, `:tensor`/`:mps`, or
`:camps`).
"""
random_product_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector) =
    product_state(rand([0,1],N); backend=backend)

"""
    zero_state(N::Int; backend=:statevector)

Return the `N`-qubit `|0…0⟩` state for the chosen backend:
`:statevector` (default, sparse vector), `:density_matrix` (sparse density
matrix), `:tensor` / `:mps` (`ITensorMPS.MPS`), or `:camps` (`CAMPSState`).
"""
function zero_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector)
    return _state_from_backend(
        backend,
        () -> sa.SparseVector(2^N, [1], [1.0+0im]),
        () -> zero_state(N_MPS(N)),
        () -> CAMPSState(N),
    )
end

"""
    one_state(N::Int; backend=:statevector)

Return the `N`-qubit `|1…1⟩` state for the chosen backend (`:statevector`,
`:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
function one_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector)
    return _state_from_backend(
        backend,
        () -> sa.SparseVector(2^N, [2^N], [1.0+0im]),
        () -> one_state(N_MPS(N)),
        () -> apply(CAMPSState(N), QuantumOps[Op("X", q) for q in 1:N]),
    )
end

"""
    plus_state(N::Int; backend=:statevector)

Return the `N`-qubit `|+…+⟩` state for the chosen backend (`:statevector`,
`:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
function plus_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector)
    return _state_from_backend(
        backend,
        () -> sa.SparseVector(fill((1.0+0im) / sqrt(2^N), 2^N)),
        () -> plus_state(N_MPS(N)),
        () -> apply(CAMPSState(N), QuantumOps[Op("H", q) for q in 1:N]),
    )
end

"""
    minus_state(N::Int; backend=:statevector)

Return the `N`-qubit `|-…-⟩` state for the chosen backend (`:statevector`,
`:density_matrix`, `:tensor`/`:mps`, or `:camps`).
"""
function minus_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector)
    return _state_from_backend(
        backend,
        () -> begin
            state = one_state(N)
            for i in 1:N
                state = Op("H", i) * state
            end
            state
        end,
        () -> minus_state(N_MPS(N)),
        () -> begin
            st = CAMPSState(N)
            for q in 1:N
                st = apply(st, Op("X", q))
                st = apply(st, Op("H", q))
            end
            st
        end,
    )
end

"""
    random_state(N::Int; backend=:statevector, chi::Int=1)

Return a random normalized `N`-qubit state for the chosen backend
(`:statevector`, `:density_matrix`, `:tensor`/`:mps`, or `:camps`). `chi` sets
the bond dimension for `:tensor` and `:camps`.
"""
function random_state(N::Int; backend::Union{Symbol,AbstractString}=:statevector, chi::Int=1)
    return _state_from_backend(
        backend,
        () -> sa.sparse(la.normalize(rand(Complex{Float64}, 2^N))),
        () -> random_state(N_MPS(N), chi),
        () -> CAMPSState(random_state(N_MPS(N), chi)),
    )
end

##========== state preparation ==========


"""
    sample(state, shots) -> Vector

Return raw computational-basis outcomes sampled from a statevector, density
matrix, or MPS. Outcomes are integers starting at 0; qubit 1 is the most
significant bit. The state is not collapsed. Use `measure(state, shots)` for
summary statistics or `measure(state)` for exact probabilities.
"""
function sample(state::AbstractVector{<:Number}, shots)

    N=get_N(state)

    # Convert state to probability distribution
    probs = abs2.(state)

    bit_basis=0:2^N-1

    sampled_outcomes = sb.sample(bit_basis, sb.Weights(probs), shots)

    # Return the sampled outcomes
    return sampled_outcomes
end

function sample(state::sa.SparseVector{<:Number}, shots)
    get_N(state) # Check the Hilbert-space dimension, as for dense states.
    indices, amplitudes = sa.findnz(state)
    # Sample only stored amplitudes, rather than scanning all 2^N basis states.
    return sb.sample(indices .- 1, sb.Weights(abs2.(amplitudes)), shots)
end


function sample_exact(state::AbstractVector{<:Number})
    probs = abs2.(state)
    indices = findall(!iszero, probs)
    return indices .- 1, probs[indices]
end

function sample_exact(state::sa.SparseVector{<:Number})
    a,b=sa.findnz(abs2.(state))
    return a .- 1,b
end
