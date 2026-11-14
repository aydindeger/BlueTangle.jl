# Statevector trajectories.

"""
`apply(state::AbstractVectorS, op::QuantumOps)`

Return the state after applying an operation. Measurements and noise sample one
trajectory. Use `run(state, op)` to return the state and measurement record.
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
`run(state, Op("MZ", qubit))` and its MX/MY variants.
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

"""
`product_state(list_of_qubits::Vector) -> sa.SparseVector`

Creates a quantum state vector from a list of qubits.

- `list_of_qubits`: A vector representing the state of each qubit.

Returns a sparse vector representing the quantum state of the system.
"""
product_state(list_of_qubits::Vector)=sa.sparse(foldl(kron,_bin2state.(Int.(sign.(list_of_qubits)))))
# product_state(N::Int,list_of_qubits::Vector)=sa.sparse(foldl(kron,_bin2state.(list_of_qubits)))
neel_state01(N::Int)=product_state([isodd(i) ? 0 : 1 for i=1:N])
neel_state10(N::Int)=product_state([isodd(i) ? 1 : 0 for i=1:N])


"""
    Z3(N::Int,j=1)
"""
Z3(N::Int,j=1)=product_state([mod(i+2-j+1,3)==0 ? 1 : 0 for i=1:N])

"""
    Z3_s(N::Int)
"""
Z3_s(N::Int)=la.normalize(Z3(N,1)+Z3(N,2)+Z3(N,3))

"""
    neel_state_s(N::Int)=(neel_state01(N)+neel_state10(N))/sqrt(2)
"""
neel_state_s(N::Int)=(neel_state01(N)+neel_state10(N))/sqrt(2)

"""
    random_product_state(N)
"""
random_product_state(N::Int)=product_state(rand([0,1],N))

"""
zero_state(N::Int) -> sa.SparseVector

Returns a sparse vector representing the |000...> quantum state.
"""
zero_state(N::Int)=sa.SparseVector(2^N, [1], [1.0+0im])
one_state(N::Int)=sa.SparseVector(2^N, [2^N], [1.0+0im])

plus_state(N::Int)=sa.SparseVector(fill((1.0+0im) / sqrt(2^N),2^N))

function minus_state(N::Int)
    state=one_state(N)
    for i=1:N
        state=Op("H",i)*state
    end
    return state
end

"""
    random_state(N::Int) -> sa.SparseVector
"""
random_state(N::Int)=sa.sparse(la.normalize(rand(Complex{Float64}, 2^N)));

##========== state preparation ==========


"""
    sample(state, shots) -> Vector

Return raw computational-basis outcomes sampled from a statevector, density
matrix, or MPS. Outcomes are integers starting at 0; qubit 1 is the most
significant bit. The state is not collapsed. Use `measure(state, shots)` for
summary statistics or `measure(state)` for exact probabilities.
"""
function sample(state::AbstractVectorS, shots)

    N=get_N(state)

    # Convert state to probability distribution
    probs = abs2.(state)

    bit_basis=0:2^N-1

    sampled_outcomes = sb.sample(bit_basis, sb.Weights(probs), shots)

    # Return the sampled outcomes
    return sampled_outcomes
end

function sample(state::sa.SparseVector, shots)
    get_N(state) # Check the Hilbert-space dimension, as for dense states.
    indices, amplitudes = sa.findnz(state)
    # Sample only stored amplitudes, rather than scanning all 2^N basis states.
    return sb.sample(indices .- 1, sb.Weights(abs2.(amplitudes)), shots)
end


function sample_exact(state::AbstractVectorS)
    probs = abs2.(state)
    indices = findall(!iszero, probs)
    return indices .- 1, probs[indices]
end

function sample_exact(state::sa.SparseVector)
    a,b=sa.findnz(abs2.(state))
    return a .- 1,b
end
