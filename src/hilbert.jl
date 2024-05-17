"""
`hilbert(N::Int, mat::AbstractMatrix, qubit::Int, target_qubit::Int)`

Constructs a sparse matrix representing the action of a quantum gate in a Hilbert space associated with a quantum system of `N` qubits

The gate `mat` is applied to the `qubit` and`target_qubit`. If `qubit` is greater than `target_qubit`, a controlled 
swap is performed before applying `mat`.

# Arguments
- `N::Int`: The number of qubits in the system.
- `mat::AbstractMatrix`: The quantum gate to be applied.
- `qubit::Int`: The qubit to which the gate is applied.
- `target_qubit::Int`: The target qubit to which the gate is applied.

# Returns
`SparseMatrix`: The resulting sparse matrix representation of the gate operation.
"""
function hilbert(N::Int,mat::AbstractMatrix,qubit::Int,target_qubit::Int;control::Int=-2)

    distance=abs(qubit-target_qubit)

    if N < qubit || N < target_qubit || N < control
        throw("N must be larger than qubits")
    elseif size(mat,1)>4
        throw("only 2-qubit operations are supported")
    end

    id = sa.sparse(sa.I, 2, 2);

    final_mat=qubit > target_qubit ? sa.sparse(BlueTangle._swap_control_target(mat)) : sa.sparse(mat)

    if control==-2

        if distance==1 #local
            e_ops = [x == qubit ? final_mat : id for x in 1:N if x!=target_qubit]
            return foldl(kron,e_ops)
        else #pauli decomposition and reconstruction to build nonlocal ops
            # println("nonlocal hilbert reconstruction")
            coefficients, _ = pauli_decomposition(final_mat)
            final_mat_nonlocal=pauli_reconstruction(coefficients,min(qubit,target_qubit);distance=distance)#nonlocal reconstruct

            for _=1:N-max(qubit,target_qubit)
                final_mat_nonlocal=kron(final_mat_nonlocal,id)
            end

            return final_mat_nonlocal

        end

    else #control operation

        if distance>1 #local
            throw("nonlocal operations (qubit and target) with control are not allowed")
        end

        list=foldl(kron,[x==qubit ? final_mat : (x==control ? sa.sparse(gate.P1) : id) for x=1:N if x!=target_qubit])    
        return list+foldl(kron,[x==control ? sa.sparse(gate.P0) : id for x=1:N])#control
    end

end



"""
`hilbert(N::Int, mat::AbstractMatrix, qubit::Int)`

Constructs a sparse matrix representing the action of a quantum gate in a Hilbert space associated with a quantum system of `N` qubits.

# Arguments
- `N::Int`: The number of qubits in the system.
- `mat::AbstractMatrix`: The quantum gate to be applied.
- `qubit::Int`: The qubit to which the gate is applied.

# Returns
`SparseMatrix`: The resulting sparse matrix representation of the gate operation.
"""
function hilbert(N::Int,mat::AbstractMatrix,qubit::Int;control::Int=-2)

    if N < qubit || N < control
        throw("N must be larger than qubit")
    elseif size(mat,1)>2 
        throw("only 2-qubit operations are supported")
    end

    id = sa.sparse(sa.I, 2, 2);

    if control==-2
        return foldl(kron,[x==qubit ? sa.sparse(mat) : id for x=1:N])
    else
        list=foldl(kron,[x==qubit ? sa.sparse(mat) : (x==control ? sa.sparse(gate.P1) : id) for x=1:N])
        return list+foldl(kron,[x==control ? sa.sparse(gate.P0) : id for x=1:N])#control
    end
end



function hilbert_layer(N::Int,layer::Vector{<:QuantumOps},state::sa.SparseVector) #todo #fix and implement control here
    
    id = sa.sparse([1.0 0; 0im 1]);
    vec=fill(id,N)

    for op=layer
        if op.q==1
            vec[op.qubit]=op.mat
        elseif op.q==2

            if abs(op.qubit-op.target_qubit)>1 #nonlocal
                throw("nonlocal operations (qubit and target) with control are not allowed")
            end

            final_mat=op.qubit > op.target_qubit ? sa.sparse(BlueTangle._swap_control_target(op.mat)) : op.mat
            vec[op.qubit]=final_mat
            vec[op.target_qubit]*=NaN
        end#
    end
    
    return foldl(kron,vec[any.(!isnan, vec)])*state

end


"""
`_swap_control_target(matrix::Matrix) -> Matrix`

Swap the control and target qubits in a two-qubit gate matrix.

- `matrix`: A 4x4 matrix representing a two-qubit operation.

Returns a new matrix with swapped control and target qubits.
"""
function _swap_control_target(matrix::Matrix)
    # Ensure the matrix is 4x4
    if size(matrix) != (4, 4)
        throw("Matrix must be 4x4")
    end

    result = zero(matrix)#zeros(Complex{Float64}, 4, 4)
    
    # Define the permutation that swaps the control and target qubits
    permutation = [1, 3, 2, 4]
    
    for i in 1:4
        for j in 1:4
            result[i, j] = matrix[permutation[i], permutation[j]]
        end
    end
    
    return result
end

# apply(op::QuantumOps,state::sa.SparseVector;noise::Union{NoiseModel,Bool}=false)=apply(state,op;noise=noise)

"""
`apply(state::sa.SparseVector, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `state`: A sparse quantum state vector to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
"""
function apply(state::sa.SparseVector,op::QuantumOps;noise::Union{NoiseModel,Bool}=false)
    
    N=get_N(state)

    # if op.q!=1 && abs(op.qubit-op.target_qubit)>1
    #     throw("non-local gate $(op.name) is not allowed!")
    # end

    if isa(op,OpF)
        state=op.apply(state)
    elseif isa(op,QC)
        state=op.apply(state)
    elseif op.type=="🔬"
        if isa(op,ifOp)
            state,ind=op.born_apply(state,noise)
        else
            state,ind=_born_measure(state,op)
        end
        # println("measurement result=$(ind)")
    else #good old gates
        state=op.expand(N)*state
    end

    ##aply noise.
    if isa(noise, NoiseModel) && op.noisy
        state=apply_noise(state,op,noise)
    end

    return state

end


"""
    apply noise on qubit or target_qubit of a given state and noise model
"""
function apply_noise(state::sa.SparseVector,op::QuantumOps,noise::NoiseModel)
    
    if op.q==1
        if op.control == -2 
            return noise.q1.apply(state,op.qubit)
        else
            return noise.q2.apply(state,op.control,op.qubit)
        end
    elseif op.q==2
        return noise.q2.apply(state,op.qubit,op.target_qubit)
    end

end



"""
apply noise on qubit or target_qubit of a given density matrix and noise model
"""
function apply_noise(rho::sa.SparseMatrixCSC,op::QuantumOps,noise::NoiseModel)
    
    if op.q==1
        if op.control == -2 
            return noise.q2.apply(rho,op.qubit)
        else
            return noise.q1.apply(rho,op.control,op.qubit)
        end
    elseif op.q==2
        return noise.q2.apply(rho,op.qubit,op.target_qubit)
    end

end




function apply(ops::Vector,psi::it.MPS;noise::Union{NoiseModel,Bool}=false,cutoff=1e-10,maxdim=500)
    for o=ops
        psi=apply(o,psi;noise=noise,cutoff=cutoff,maxdim=maxdim)
    end

    return psi
end

apply(op::QuantumOps,state::it.MPS;noise::Union{NoiseModel,Bool}=false,cutoff=1e-10,maxdim=500)=apply(state,op;noise=noise,cutoff=cutoff,maxdim=maxdim)

function apply(psi::it.MPS,op::QuantumOps;noise::Union{NoiseModel,Bool}=false,cutoff=1e-10,maxdim=500)
    
    M=get_N(psi)

    if op.q!=1 && abs(op.qubit-op.target_qubit)>1
        throw("non-local gate $(op.name) is not allowed!")
    end

    if op.type=="🔬"
        if isa(op,ifOp)
            throw("MPS ifOp measurements does not work yet.")
            # state,ind=op.born_apply(state,noise) #fix
        else
            psi,ind=_born_measure(psi,op;cutoff=cutoff,maxdim=maxdim)
        end
    else
        psi=it.apply(op.expand(M),psi;cutoff=cutoff,maxdim=maxdim)
    end

    if isa(noise, NoiseModel) && op.noisy
        throw("Noisy MPS does not work yet.")
    #     selected_noise = op.q == 1 ? noise.q1 : noise.q2
    #     if isa(selected_noise, QuantumChannel)
    #         psi = apply_noise(psi, op, selected_noise)
    #     end
    end
    
    return psi

end


"""
`apply(rho::sa.SparseMatrixCSC, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `rho`: A sparse quantum density matrix to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
"""
function apply(rho::sa.SparseMatrixCSC,op::QuantumOps;noise::Union{NoiseModel,Bool}=false)

    N=get_N(rho)
    
    if op.q!=1 && abs(op.qubit-op.target_qubit)>1
        throw("non-local gate $(op.name) is not allowed! Use different layout or compile circuit to add swaps")
    end

    if op.type=="🔬"
        isa(op,ifOp) ? rho=op.born_apply(rho,noise) : rho=_born_measure(rho,op)
    else
        e_op=op.expand(N)
        rho=e_op*rho*e_op'
    end

    if isa(noise, NoiseModel) && op.noisy
        selected_noise = op.q == 1 ? noise.q1 : noise.q2
        if isa(selected_noise, QuantumChannel)
            rho = apply_noise(rho, op, selected_noise)
        end
    end

    return rho
    
end

function _born_measure(state::sa.SparseVector,o::QuantumOps)

    N=get_N(state)
    rotMat=o.expand(N)
    state=rotMat*state#rotate
    state,ind=born_measure_Z(N,state,o.qubit)
    state=rotMat'*state#rotate back

    return state,ind

end


function born_measure_Z(N::Int,state::sa.SparseVector,qubit::Int)

    born_ops=[gate.P0, gate.P1]

    #according to benchmark
    if N<=12
        prob0=sum(abs2.(hilbert(N,born_ops[1],qubit)*state))
    else
        prob0=real(BlueTangle.partial_trace(state,qubit))[1,1]
    end

    ind=rand() < prob0 ? 0 : 1
    return sa.normalize(hilbert(N,born_ops[ind+1],qubit)*state),ind

end


"""
    born_measure_Z(N::Int,state::sa.SparseVector,qubit1::Int,qubit2::Int)
    born_measure_Z(N::Int,state::sa.SparseVector,qubit1::Int)

    #todo test this!
"""
function born_measure_Z(N::Int,state::sa.SparseVector,qubit1::Int,qubit2::Int)
  
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



function _born_measure(psi::it.MPS,o::QuantumOps;cutoff=1e-10,maxdim=500)

    M=it.siteinds(psi)
    psi=it.apply(o.expand(M),psi;cutoff=cutoff,maxdim=maxdim) #rotate
    psi,ind=born_measure_Z(psi,o.qubit)
    psi=it.apply(it.op(o.mat',M[o.qubit]),psi;cutoff=cutoff,maxdim=maxdim)#rotate back
    return psi,ind

end

function born_measure_Z(psi::it.MPS,qubit::Int;cutoff=1e-10,maxdim=500)

    born_ops=[gate.P0, gate.P1]
    si=it.siteinds(psi, qubit)

    it.orthogonalize!(psi, qubit)
    psij=psi[qubit]
    rho_j = it.prime(it.dag(psij), si) * psij
    prob0=real.(rho_j[1])#0
    ind=rand() < prob0 ? 0 : 1

    psi=it.normalize(it.apply(it.op(born_ops[ind+1],si),psi;cutoff=cutoff,maxdim=maxdim))

    return psi,ind

end

function _born_measure(rho::sa.SparseMatrixCSC,o::QuantumOps)

    throw("fix this:")
    N=get_N(rho)
    
    mat = o.expand(N)
    rho = mat * rho * mat'; #rotate

    rho_born =born_measure_Z(N,rho,o.qubit)
    
    return mat' * rho_born * mat#todo fix rotate back

end

function born_measure_Z(N::Int,rho::sa.SparseMatrixCSC,qubit::Int)

    new_rho=zero(rho)

    born_ops=[gate.P0, gate.P1]
    for b=born_ops
        eb=hilbert(N,b,qubit)
        new_rho += eb * rho * eb'
    end
    
    return new_rho

end

"""
`_weighted_sample(ek_ops::Vector{Any}, probs::Vector{Float64}) -> Tuple`

Selects a Kraus operator from a set based on their associated probabilities.

- `ek_ops`: Vector of extended Kraus operators.
- `probs`: Vector of probabilities corresponding to each Kraus operator.

Returns a tuple containing the index and the selected Kraus operator.
"""

# function _weighted_sample(probs::Vector{Any})
function _weighted_sample(probs::Vector)

    rval=rand()

    for (i, cum_weight) in enumerate(cumsum(probs))
        if rval <= cum_weight
            return i
        end
    end
end



##
##========== state preparation ==========

"""
`product_state(list_of_qubits::Vector) -> sa.SparseVector`

Creates a quantum state vector from a list of qubits.

- `list_of_qubits`: A vector representing the state of each qubit.

Returns a sparse vector representing the quantum state of the system.
"""
product_state(list_of_qubits::Vector)=sa.sparse(foldl(kron,_bin2state.(list_of_qubits)))
product_state(N::Int,list_of_qubits::Vector)=sa.sparse(foldl(kron,_bin2state.(list_of_qubits)))
neel_state01(N::Int)=product_state([isodd(i) ? 0 : 1 for i=1:N])
neel_state10(N::Int)=product_state([isodd(i) ? 1 : 0 for i=1:N])

"""
    random_state(N)
"""
random_state(N::Int)=product_state(rand([0,1],N))

"""
zero_state(N::Int) -> sa.SparseVector

Returns a sparse vector representing the |000...> quantum state.
"""
zero_state(N::Int)=sa.SparseVector(2^N, [1], [1.0+0im])
one_state(N::Int)=sa.SparseVector(2^N, [2^N], [1.0+0im])

"""
    mixed_state(N::Int) -> sa.SparseVector
"""
mixed_state(N::Int)=sa.sparse(la.normalize(rand(Complex{Float64}, 2^N)));

##========== state preparation ==========