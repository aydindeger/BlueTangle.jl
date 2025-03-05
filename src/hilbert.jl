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
    # println("N=$(N), qubit=$(qubit), target_qubit=$(target_qubit)")
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

    else #control operation on two-qubit gates
            
        if distance==1 #nonlocal #check this
            list=foldl(kron,[x==qubit ? final_mat : (x==control ? sa.sparse(gate.P1) : id) for x=1:N if x!=target_qubit])    
            return list+foldl(kron,[x==control ? sa.sparse(gate.P0) : id for x=1:N])#control
            
        else

            if mat==gate.CX
                return CCZX("CCX",N,qubit,control,target_qubit)
            elseif mat==gate.CZ
                return CCZX("CCZ",N,qubit,control,target_qubit)
            else
                throw("Unsupported operation. Use 'CCZ' or 'CCX'.")
            end
            
        end

    end

end


function CCZX(sym::String,N::Int,i::Int,j::Int,k::Int)

    if any(q -> q < 1 || q > N, [i, j, k])
        throw("Qubit indices must be within the range 1 to N")
    end

    i=N-i
    j=N-j
    k=N-k

    H = sa.spzeros(2^N,2^N);

    if sym=="CCZ"
        for a in 0:2^N-1
            H[a+1,a+1]+= ((a >> i) & (a >> j) & (a >> k) & 1)==1 ? -1 : 1 
        end
    elseif sym=="CCX"
        for a in 0:2^N-1
            if ((a >> i) & (a >> j) & 1)==1
                b=a⊻(1<<k)
                H[a+1,b+1]+=1;
            else
                H[a+1,a+1]+=1;
            end
        end
    else
        throw("Unsupported operation. Use 'CCZ' or 'CCX'.")
    end
    
    return H
end


function hilbert3(N::Int,mat::AbstractMatrix,first_qubit::Int)

    second_qubit=first_qubit+1
    third_qubit=second_qubit+1

    if N < third_qubit
        throw("N must be larger than all three qubits")
    elseif size(mat,1)!=8
        throw("only 3-qubit operations are supported")
    end

    id = sa.sparse(ComplexF64,sa.I, 2, 2);
    e_ops=fill(id,N)

    insert!(e_ops,first_qubit,mat)

    for i=1:3
        popat!(e_ops,first_qubit+1)
    end

    return foldl(kron,e_ops)

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



# function hilbert_layer(N::Int,layer::Vector{<:QuantumOps},state::AbstractVectorS) #todo #fix and implement control here
    
#     id = sa.sparse([1.0 0; 0im 1]);
#     vec=fill(id,N)

#     for op=layer
#         if op.q==1
#             vec[op.qubit]=op.mat
#         elseif op.q==2

#             if abs(op.qubit-op.target_qubit)>1 #nonlocal
#                 throw("nonlocal operations (qubit and target) with control are not allowed")
#             end

#             final_mat=op.qubit > op.target_qubit ? sa.sparse(BlueTangle._swap_control_target(op.mat)) : op.mat
#             vec[op.qubit]=final_mat
#             vec[op.target_qubit]*=NaN
#         end#
#     end
    
#     return foldl(kron,vec[any.(!isnan, vec)])*state

# end


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

"""
    relabel_swap(ops::Vector{Op}) -> Tuple{Vector{Op}, Vector{Int}}

Relabel qubit indices in a sequence of quantum operations, simulating the effect of SWAP gates without actually performing them.

This function processes a list of quantum operations, updating qubit indices to reflect the effect of SWAP operations. It handles both SWAP gates and other quantum operations, adjusting their qubit references according to the cumulative effect of all SWAP operations.

# Arguments
- `ops::Vector{Op}`: A vector of `Op` objects representing quantum operations.

# Returns
- `Tuple{Vector{Op}, Vector{Int}}`: A tuple containing:
  1. A vector of relabeled `Op` objects, excluding SWAP operations.
  2. The final qubit mapping, where `mapping[i]` gives the final label of the qubit initially labeled `i`.

# Details
- SWAP operations are used to update the qubit mapping but are not included in the output operations.
- For non-SWAP operations, qubit indices are updated based on the current mapping.
- The function assumes qubit indices start from 1.

relabeled_ops, final_mapping = relabel_swap(ops)
```
"""
function relabel_swap(ops::Vector{T}) where T<:QuantumOps
    N = maximum(max(op.qubit, op.target_qubit, op.control) for op in ops)
    qubit_mapping = collect(1:N)
    relabeled_ops = Op[]

    for op in ops
        if op.name == "SWAP"

            if op.control!=-2
                throw("Swap relabeling doesn't work with the control qubit.")
            end

            i, j = op.qubit, op.target_qubit
            # Relabel the SWAP operation
            new_i = findfirst(==(i), qubit_mapping)
            new_j = findfirst(==(j), qubit_mapping)

            qubit_mapping[new_i], qubit_mapping[new_j] = qubit_mapping[new_j], qubit_mapping[new_i]

        else
            # Relabel the operation based on current mapping
            new_qubit = findfirst(==(op.qubit), qubit_mapping)
            new_target = op.target_qubit > 0 ? findfirst(==(op.target_qubit), qubit_mapping) : op.target_qubit
            new_control = op.control > 0 ? findfirst(==(op.control), qubit_mapping) : op.control

            new_op=Op(op.name,new_qubit,new_target;control=new_control,noisy=op.noisy)

            push!(relabeled_ops, new_op)
        end
    end

    return relabeled_ops, qubit_mapping
end

"""
    apply noise on qubit or target_qubit of a given state and noise model
"""
function apply_noise(state::AbstractVectorS,op::QuantumOps,noise::NoiseModel)
    
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
            return noise.q1.apply(rho,op.qubit)
        else
            return noise.q2.apply(rho,op.control,op.qubit)
        end
    elseif op.q==2
        return noise.q2.apply(rho,op.qubit,op.target_qubit)
    end

end

# apply(op::QuantumOps,state::AbstractVectorS;noise::Union{NoiseModel,Bool}=false)=apply(state,op;noise=noise)

"""
`apply(state::AbstractVectorS, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `state`: A sparse quantum state vector to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
"""
function apply(state::AbstractVectorS,op::QuantumOps;noise::Union{NoiseModel,Bool}=false)
    
    N=get_N(state)

    # if op.q!=1 && abs(op.qubit-op.target_qubit)>1
    #     throw("non-local gate $(op.name) is not allowed!")
    # end

    if isa(op,OpF)
        state=op.apply(state)
    elseif isa(op,OpQC)
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

function apply(psi::it.MPS,op::QuantumOps;noise::Union{NoiseModel,Bool}=false,cutoff=1e-12)#,cutoff=1e-10,maxdim=500)
    
    M=get_M(psi)

    # if op.q!=1 && abs(op.qubit-op.target_qubit)>1
    #     throw("non-local gate $(op.name) is not supported.")
    # end

    if isa(op,OpF)
        psi=op.apply(psi;cutoff=cutoff)
    elseif isa(op,OpQC)
        throw("Quantum Channel MPS is not supported")
        # psi=op.apply(psi)
    elseif op.type=="🔬"
        if isa(op,ifOp)
            throw("MPS ifOp measurement is not supported")
            # psi,ind=op.born_apply(psi,noise)
        else
            psi,ind=_born_measure(psi,op)#;cutoff=cutoff,maxdim=maxdim)
        end
        # println("measurement result=$(ind)")
    else #good old gates
        psi=it.normalize(it.apply(op.expand(M),psi;cutoff=cutoff))#;cutoff=cutoff,maxdim=maxdim))
    end

    ##aply noise.
    if isa(noise, NoiseModel) && op.noisy
        throw("Noisy MPS is not supported")
        # state=apply_noise(state,op,noise)
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

    if isa(op,OpF)
        rho=op.apply(rho)
    elseif isa(op,OpQC)
        rho=op.apply(rho)
    elseif op.type=="🔬"
        if isa(op,ifOp)
            rho,ind=op.born_apply(rho,noise)
        else
            rho,ind=_born_measure(rho,op)
        end
        # println("measurement result=$(ind)")
    else #good old gates
        e_op=op.expand(N)
        rho=e_op*rho*e_op'
    end

    ##aply noise
    if isa(noise, NoiseModel) && op.noisy
        rho=apply_noise(rho,op,noise)
    end

    return rho
    
end


function _born_measure(state::AbstractVectorS,o::QuantumOps)

    N=get_N(state)
    rotMat=o.expand(N)
    state=rotMat*state#rotate
    state,ind=born_measure_Z(N,state,o.qubit)
    state=rotMat'*state#rotate back

    return state,ind

end


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
    born_measure_Z(N::Int,state::AbstractVectorS,qubit1::Int)

    #todo test this!
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



function _born_measure(psi::it.MPS,o::QuantumOps)#;cutoff=1e-10,maxdim=500)

    M=it.siteinds(psi)
    psi=it.normalize(it.apply(o.expand(M),psi))#;cutoff=cutoff,maxdim=maxdim)) #rotate
    psi,ind=born_measure_Z(psi,o.qubit)
    psi=it.normalize(it.apply(it.op(o.mat',M[o.qubit]),psi))#;cutoff=cutoff,maxdim=maxdim))#rotate back
    return psi,ind

end

function born_measure_Z(psi::it.MPS,qubit::Int)#;cutoff=1e-10,maxdim=500)

    born_ops=[gate.P0, gate.P1]
    si=it.siteinds(psi, qubit)

    it.orthogonalize!(psi, qubit)
    psij=psi[qubit]
    rho_j = it.prime(it.dag(psij), si) * psij
    prob0=real.(rho_j[1])#0
    ind=rand() < prob0 ? 0 : 1

    psi=it.normalize(it.apply(it.op(born_ops[ind+1],si),psi))#;cutoff=cutoff,maxdim=maxdim))

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