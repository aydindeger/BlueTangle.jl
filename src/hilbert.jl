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
function hilbert(N::Int, mat::AbstractMatrix, qubit::Int, target_qubit::Int; control::Int=-2)
    all(q -> 1 <= q <= N, (qubit, target_qubit)) ||
        throw(ArgumentError("qubit indices must be within 1:N"))
    qubit != target_qubit || throw(ArgumentError("operated qubits must differ"))
    size(mat) == (4, 4) || throw(ArgumentError("a two-qubit matrix must be 4×4"))
    if control != -2
        1 <= control <= N && control != qubit && control != target_qubit ||
            throw(ArgumentError("control must be a distinct qubit within 1:N"))
        # The control projector commutes with the two-qubit operator because
        # their supports are disjoint. This also handles nonlocal custom gates.
        return hilbert(N, gate.P0, control) +
               hilbert(N, gate.P1, control) * hilbert(N, mat, qubit, target_qubit)
    end

    distance = abs(qubit - target_qubit)
    identity = sa.sparse(la.I, 2, 2)
    ordered_mat = qubit > target_qubit ? sa.sparse(_swap_control_target(mat)) : sa.sparse(mat)
    if distance == 1
        return foldl(kron, [q == qubit ? ordered_mat : identity for q in 1:N if q != target_qubit])
    end
    coefficients, _, _ = pauli_decomposition(ordered_mat)
    expanded = pauli_reconstruction(coefficients, min(qubit, target_qubit); distance=distance)
    for _ in 1:N-max(qubit, target_qubit)
        expanded = kron(expanded, identity)
    end
    return expanded
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




"""
the gate is invariant under swapping the qubits.
"""
function _check_swap_invariant(matrix::Union{AbstractMatrixS,Function})

    if isa(matrix,Function)
        arg_no=BlueTangle._find_argument_number(matrix)
        matrix=matrix(la.rand(arg_no)...)
    end

    if size(matrix) != (4, 4)
        throw("Matrix must be 4x4")
    end

    return all(x->x==0,matrix*gate.SWAP-gate.SWAP*matrix)

end

"""
`_swap_control_target(matrix::Matrix) -> Matrix`

Swap the control and target qubits in a two-qubit gate matrix.

- `matrix`: A 4x4 matrix representing a two-qubit operation.

Returns a new matrix with swapped control and target qubits.
"""
function _swap_control_target(matrix::Union{AbstractMatrixS,Function})
    
    if _check_swap_invariant(matrix)
        return matrix
    end

    if isa(matrix,Function)
        throw("this function is not invariant under swap and it's a problem!")
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
    relabel_swap(ops::Vector{Op}, qubit_mapping::Vector{Int}=Int[]) -> Tuple{Vector{Op}, Vector{Int}}

Relabel qubit indices in a sequence of quantum operations, simulating the effect of SWAP gates without actually performing them.

This function processes a list of quantum operations, updating qubit indices to reflect the effect of SWAP operations. It handles both SWAP gates and other quantum operations, adjusting their qubit references according to the cumulative effect of all SWAP operations.

# Arguments
- `ops::Vector{Op}`: A vector of `Op` objects representing quantum operations.
- `qubit_mapping::Vector{Int}=Int[]`: Existing qubit mapping to continue from. If omitted,
  the mapping starts from the identity permutation.

# Returns
- `Tuple{Vector{Op}, Vector{Int}}`: A tuple containing:
  1. A vector of relabeled `Op` objects, excluding SWAP operations.
  2. The final qubit mapping, where `mapping[i]` gives the final label of the qubit initially labeled `i`.

# Details
- SWAP operations are used to update the qubit mapping but are not included in the output operations.
- For non-SWAP operations, qubit indices are updated based on the current mapping.
- The function assumes qubit indices start from 1 and that any provided `qubit_mapping`
  is a permutation of `1:N`. If later operations reference higher-index qubits that were
  not tracked yet, the mapping is extended with identity entries for those qubits.

relabeled_ops, final_mapping = relabel_swap(ops)
```
"""
function _relabel_operation(op::QuantumOps, mapping::Vector{Int})
    relabel(q) = q > 0 ? findfirst(==(q), mapping) : q
    if op isa Op
        return Op(op.name, op.mat, relabel(op.qubit), relabel(op.target_qubit);
                  control=relabel(op.control), noisy=op.noisy)
    elseif op isa OpQC
        if op.q == 3
            qubits = relabel.(op.qubit:op.qubit+2)
            qubits == collect(first(qubits):first(qubits)+2) ||
                throw(ArgumentError("three-qubit channels require consecutive qubits in their original order after relabeling"))
        end
        return OpQC(op.name, op.kraus, relabel(op.qubit), relabel(op.target_qubit); type=op.type)
    elseif op isa ifOp
        return ifOp(op.name, relabel(op.qubit),
                    Op[_relabel_operation(o, mapping) for o in op.if01[1]],
                    Op[_relabel_operation(o, mapping) for o in op.if01[2]])
    elseif op isa OpF && op.data isa Vector{<:QuantumOps}
        return OpF(op.name, [_relabel_operation(o, mapping) for o in op.data])
    elseif mapping == collect(eachindex(mapping))
        return op
    end
    throw(ArgumentError("an opaque OpF cannot be relabeled; provide its constituent gates instead"))
end

function relabel_swap(ops::Vector{T}, qubit_mapping::Vector{Int}=Int[]) where T<:QuantumOps
    max_qubit = get_N(ops)
    if isempty(qubit_mapping)
        qubit_mapping = collect(1:max_qubit)
    else
        N = length(qubit_mapping)
        sort(qubit_mapping) == collect(1:N) ||
            throw(ArgumentError("qubit_mapping must be a permutation of 1:$N"))
        qubit_mapping = copy(qubit_mapping)
        max_qubit > N && append!(qubit_mapping, (N + 1):max_qubit)
    end

    relabeled_ops = T[]
    for op in ops
        # Only an actual, unconditional SWAP matrix may be absorbed into labels.
        if op isa Op && uppercase(op.name) == "SWAP" && op.mat == gate.SWAP
            op.control == -2 ||
                throw(ArgumentError("a controlled SWAP cannot be replaced by a fixed relabeling"))
            i = findfirst(==(op.qubit), qubit_mapping)
            j = findfirst(==(op.target_qubit), qubit_mapping)
            qubit_mapping[i], qubit_mapping[j] = qubit_mapping[j], qubit_mapping[i]
        else
            push!(relabeled_ops, _relabel_operation(op, qubit_mapping))
        end
    end
    return relabeled_ops, qubit_mapping
end

"""
    apply noise on qubit or target_qubit of a given state and noise model
"""
function apply_noise(state::Union{AbstractVectorS,AbstractMatrixS}, op::QuantumOps, noise::NoiseModel)
    hasproperty(op, :noisy) && op.noisy || return state
    if op.q == 1 && op.control == -2
        return noise.q1 isa QuantumChannel ? noise.q1.apply(state, op.qubit) : state
    elseif op.q == 1
        return noise.q2 isa QuantumChannel ? noise.q2.apply(state, op.control, op.qubit) : state
    elseif op.q == 2
        return noise.q2 isa QuantumChannel ? noise.q2.apply(state, op.qubit, op.target_qubit) : state
    end
    return state
end


function _apply_quantum_channel_mps(psi::it.MPS, channel::QuantumChannel, qubits::Int...; kwargs...)
    return _apply_quantum_channel_mps(psi, channel.kraus, qubits...; kwargs...)
end

function _apply_quantum_channel_mps(psi::it.MPS, kraus_ops::Vector{<:AbstractMatrix},
                                    qubits::Int...; kwargs...)
    sites = get_M(psi)
    # ITensor orders local matrix indices from least to most significant.
    local_sites = sites[reverse(collect(qubits))]
    tensors = [it.op(K, local_sites...) for K in kraus_ops]
    # Evaluate Born probabilities without truncation, and without expanding the
    # MPS into a 2^N state vector. User truncation applies to the selected branch.
    probabilities = [la.norm(it.apply(K, psi; cutoff=0.0))^2 for K in tensors]
    index = _weighted_sample(probabilities)
    return la.normalize(it.apply(tensors[index], psi; kwargs...))
end


function apply_noise(psi::it.MPS,op::QuantumOps,noise::NoiseModel;kwargs...)

    if !(hasproperty(op,:noisy) && op.noisy==true)
        return psi
    end

    if op.q==1
        if op.control == -2
            if isa(noise.q1,QuantumChannel)
                return _apply_quantum_channel_mps(psi,noise.q1,op.qubit;kwargs...)
            end
        else
            if isa(noise.q2,QuantumChannel)
                return _apply_quantum_channel_mps(psi,noise.q2,op.control,op.qubit;kwargs...)
            end
        end
    elseif op.q==2
        if isa(noise.q2,QuantumChannel)
            return _apply_quantum_channel_mps(psi,noise.q2,op.qubit,op.target_qubit;kwargs...)
        end
    end

    return psi
end



# apply(op::QuantumOps,state::AbstractVectorS;noise::Union{NoiseModel,Bool}=false)=apply(state,op;noise=noise)

"""
`apply(state::AbstractVectorS, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `state`: A sparse quantum state vector to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
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
                state, mid_measurements = apply(op.data, state; noise=noise, track_measurements=true, kwargs...)
            else
                state = apply(op.data, state; noise=noise, kwargs...)
            end
        else
            state = op.apply(state)
        end
    elseif isa(op,OpQC)
        if uppercase(op.name)=="RES" || uppercase(op.name)=="RESET"
            state,_=_reset_Z(state,op.qubit)
        else
            state=op.apply(state)
        end
    elseif op.type=="🔬"
        if isa(op,ifOp)
            state,ind=op.born_apply(state,noise)
        else
            state,ind=_born_measure(state,op)
        end
        track_measurements && push!(mid_measurements,ind)
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

function apply(
    vector_of_ops::Union{Vector{<:QuantumOps},AbstractVector{<:Tuple}},
    psi::Union{AbstractVectorS,AbstractMatrixS,it.MPS};
    noise::Union{NoiseModel,Bool}=false,
    track_measurements::Bool=false,
    kwargs...
)
    mid_measurements=Int[]
    
    if isa(vector_of_ops, Vector{<:QuantumOps})

        for op in vector_of_ops
            if track_measurements
                psi,new_vals = apply(psi, op; noise=noise, track_measurements=true, kwargs...)
                append!(mid_measurements,new_vals)
            else
                psi = apply(psi, op; noise=noise, kwargs...)
            end
        end

    elseif isa(vector_of_ops, AbstractVector{<:Tuple})

        for str in vector_of_ops
            if track_measurements
                psi,new_vals = apply(psi, Op(str...); noise=noise, track_measurements=true, kwargs...)
                append!(mid_measurements,new_vals)
            else
                psi = apply(psi, Op(str...); noise=noise, kwargs...)
            end
        end

    else
        throw(ArgumentError("ops are not valid!"))
    end

    return track_measurements ? (psi,mid_measurements) : psi
end

apply(op::QuantumOps, psi::Union{it.MPS,AbstractVectorS,AbstractMatrixS}; noise::Union{NoiseModel,Bool}=false, kwargs...) = apply(psi, op; noise=noise, kwargs...)

apply(str::Tuple, psi::Union{it.MPS,AbstractVectorS,AbstractMatrixS}; noise::Union{NoiseModel,Bool}=false, kwargs...) = apply(psi, Op(str...); noise=noise, kwargs...)

function _unsupported_keyword_check(kwargs, keys_to_check::Vector)
    invalid_keys = filter(k -> k ∉ keys_to_check, keys(kwargs))
    if !isempty(invalid_keys)
        throw("Unsupported keyword argument(s) provided: $(invalid_keys)")
    end
end

function apply(
    psi::it.MPS,
    op::QuantumOps;
    noise::Union{NoiseModel,Bool}=false,
    track_measurements::Bool=false,
    kwargs...
)#,cutoff=1e-10,maxdim=500)
    
    _unsupported_keyword_check(kwargs,[:cutoff,:maxdim,:pars])

    M=get_M(psi)
    mid_measurements=Int[]
    # Decide which parameter to pass: if cutoff is defined then use it, otherwise use maxdim.
    dim_kw = Dict{Symbol,Any}()
    if haskey(kwargs, :cutoff)
        dim_kw[:cutoff] = kwargs[:cutoff]
    end
    if haskey(kwargs, :maxdim)
        dim_kw[:maxdim] = kwargs[:maxdim]
    end
    dim_kw = (; dim_kw...) # Convert dictionary to a named tuple for kwarg splatting

    # if op.q!=1 && abs(op.qubit-op.target_qubit)>1
    #     throw("non-local gate $(op.name) is not supported.")
    # end

    if isa(op,OpF)
        if op.data isa Vector{<:QuantumOps}
            if track_measurements
                psi, mid_measurements = apply(op.data, psi; noise=noise, track_measurements=true, dim_kw...)
            else
                psi = apply(op.data, psi; noise=noise, dim_kw...)
            end
        else
            psi = op.apply(psi; dim_kw...)
        end
    elseif isa(op,OpQC)
        if uppercase(op.name)=="RES" || uppercase(op.name)=="RESET"
            psi,_ = _reset_Z(psi,op.qubit)
        elseif op.q==1
            psi = _apply_quantum_channel_mps(psi,op.kraus,op.qubit; dim_kw...)
        elseif op.q==2
            psi = _apply_quantum_channel_mps(psi,op.kraus,op.qubit,op.target_qubit; dim_kw...)
        elseif op.q==3
            psi = _apply_quantum_channel_mps(psi,op.kraus,op.qubit,op.qubit+1,op.qubit+2; dim_kw...)
        end
    elseif op.type=="🔬"
        if isa(op,ifOp)
            psi,ind=op.born_apply(psi,noise)
        else
            psi,ind=_born_measure(psi,op)
        end
        track_measurements && push!(mid_measurements,ind)
        # println("measurement result=$(ind)")
    elseif op.type=="f" #for variational circuits
        psi=la.normalize(it.apply(op.expand(M,kwargs[:pars]...),psi; dim_kw...))
    else #good old gates
        psi=la.normalize(it.apply(op.expand(M),psi; dim_kw...))
    end

    ##apply noise.
    if isa(noise, NoiseModel)
        psi=apply_noise(psi,op,noise; dim_kw...)
    end
    
    return track_measurements ? (psi,mid_measurements) : psi

end



"""
    apply(rho::AbstractMatrix, op::QuantumOps; noise=false, track_measurements=false)

Evolve a density matrix. By default, measurements and conditional operations
sum their unnormalized branches, and `MR` averages over X, Y and Z bases.
With `track_measurements=true`, sample measurement outcomes instead and return
`(conditional_density_matrix, recorded_bits)`. Noise channels always apply the
complete Kraus sum. The input matrix is not modified.
"""
function apply(rho::AbstractMatrixS, op::QuantumOps;
               noise::Union{NoiseModel,Bool}=false, track_measurements::Bool=false, kwargs...)
    _unsupported_keyword_check(kwargs, [:pars])
    N = get_N(rho)
    bits = Int[]
    if op isa OpF
        if op.data isa Vector{<:QuantumOps}
            if track_measurements
                rho, bits = apply(op.data, rho; noise=noise, track_measurements=true, kwargs...)
            else
                rho = apply(op.data, rho; noise=noise, kwargs...)
            end
        else
            rho = op.apply(rho)
        end
    elseif op isa OpQC
        rho = op.apply(rho)
    elseif op.type == "🔬"
        if op isa ifOp
            rho, bit = op.born_apply(rho, noise; selective=track_measurements)
        else
            rho, bit = _measure_density(rho, op.name, op.qubit; selective=track_measurements)
        end
        track_measurements && push!(bits, bit)
    else
        operator = op.type == "f" ? op.expand(N, kwargs[:pars]...) : op.expand(N)
        rho = operator * rho * operator'
    end
    noise isa NoiseModel && (rho = apply_noise(rho, op, noise))
    return track_measurements ? (rho, bits) : rho
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

    The second return value is a 1-based outcome index ordered as 00, 10, 01, 11.
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



function _born_measure(psi::it.MPS, o::QuantumOps)
    sites = get_M(psi)
    rotation = _measurement_mat(_resolve_measurement_name(o.name))
    psi = la.normalize(it.apply(it.op(rotation, sites[o.qubit]), psi))
    psi, bit = born_measure_Z(psi, o.qubit)
    psi = la.normalize(it.apply(it.op(rotation', sites[o.qubit]), psi))
    return psi, bit
end


function born_measure_Z(psi::it.MPS,qubit::Int)#;cutoff=1e-10,maxdim=500)

    born_ops=[gate.P0, gate.P1]
    si=get_M(psi)[qubit]

    it.orthogonalize!(psi, qubit)
    psij=psi[qubit]
    rho_j = it.prime(it.dag(psij), si) * psij
    prob0=real.(rho_j[1])#0
    ind=rand() < prob0 ? 0 : 1

    psi=la.normalize(it.apply(it.op(born_ops[ind+1],si),psi))#;cutoff=cutoff,maxdim=maxdim))

    return psi,ind

end

function _reset_Z(state::AbstractVectorS,qubit::Int)
    N=get_N(state)
    state,ind=born_measure_Z(N,state,qubit)
    if ind==1
        state=Op("X",qubit)*state
    end
    return state,ind
end

function _reset_Z(psi::it.MPS,qubit::Int)
    psi,ind=born_measure_Z(psi,qubit)
    if ind==1
        M=get_M(psi)
        psi=la.normalize(it.apply(it.op(gate.X,M[qubit]),psi))
    end
    return psi,ind
end

# Project in the physical measurement basis. Leaving branches unnormalized
# keeps both Born probabilities available for deterministic ensemble evolution.
function _density_measurement_branches(rho::AbstractMatrixS, name::String, qubit::Int)
    rotation = _measurement_mat(name)
    N = get_N(rho)
    return map((gate.P0, gate.P1)) do projector
        extended = hilbert(N, rotation' * projector * rotation, qubit)
        extended * rho * extended'
    end
end

function _measure_density(rho::AbstractMatrixS, name::String, qubit::Int; selective=false)
    if selective
        branches = _density_measurement_branches(rho, _resolve_measurement_name(name), qubit)
        probabilities = max.(0.0, real.(la.tr.(branches)))
        total = sum(probabilities)
        total > 0 || throw(ArgumentError("cannot sample a density matrix with zero trace"))
        bit = rand() < probabilities[1] / total ? 0 : 1
        return branches[bit + 1] / probabilities[bit + 1], bit
    end

    bases = uppercase(name) in ("MR", "M(R)") ? ("MX", "MY", "MZ") : (name,)
    result = zero(rho)
    for basis in bases
        result += sum(_density_measurement_branches(rho, basis, qubit)) / length(bases)
    end
    return result, nothing
end

_born_measure(rho::AbstractMatrixS, op::QuantumOps) = first(_measure_density(rho, op.name, op.qubit))

"""
    born_measure_Z(N, rho::AbstractMatrix, qubit)

Return the unconditioned Z-measurement channel, retaining both outcomes.
"""
born_measure_Z(N::Int, rho::AbstractMatrixS, qubit::Int) =
    first(_measure_density(rho, "MZ", qubit))


"""
`_weighted_sample(ek_ops::Vector{Any}, probs::Vector{Float64}) -> Tuple`

Selects a Kraus operator from a set based on their associated probabilities.

- `ek_ops`: Vector of extended Kraus operators.
- `probs`: Vector of probabilities corresponding to each Kraus operator.

Returns a tuple containing the index and the selected Kraus operator.
"""

function _weighted_sample(probs::Vector)
    all(isfinite, probs) || throw(ArgumentError("probabilities must be finite"))
    all(p -> p >= -1e-12, probs) || throw(ArgumentError("probabilities must be nonnegative"))
    weights = max.(0.0, probs)
    total = sum(weights)
    total > 0 || throw(ArgumentError("at least one probability must be positive"))
    threshold = rand() * total
    cumulative = 0.0
    for i in eachindex(weights)
        cumulative += weights[i]
        threshold < cumulative && return i
    end
    return findlast(>(0), weights) # Protect against rounding at the upper endpoint.
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
