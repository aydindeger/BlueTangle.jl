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



"""
    apply(state, ops; noise=false, track_measurements=false, kwargs...)

Apply operations in their listed order to a statevector, density matrix, or MPS.
No routing or scheduling is performed. With `track_measurements=true`, return
the final state and the recorded mid-circuit outcomes.
"""
function apply(
    psi::Union{AbstractVector{<:Number},AbstractMatrixS,it.MPS},
    vector_of_ops::Union{Vector{<:QuantumOps},AbstractVector{<:Tuple}};
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

function _unsupported_keyword_check(kwargs, keys_to_check::Vector)
    invalid_keys = filter(k -> k ∉ keys_to_check, keys(kwargs))
    if !isempty(invalid_keys)
        throw("Unsupported keyword argument(s) provided: $(invalid_keys)")
    end
end

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
