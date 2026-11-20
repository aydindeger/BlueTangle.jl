"""
    hilbert(N::Int, mat::AbstractMatrix, qubit::Int, target_qubit::Int; control=-2)

Embed a 4×4 matrix on two distinct qubits of an `N`-qubit register, returning
a sparse matrix. The matrix's first tensor factor acts on `qubit` and its
second on `target_qubit`, including when their indices are reversed or nonlocal.
Qubit 1 is the most significant tensor factor of the full register.

Set `control` to a third, distinct qubit to apply `mat` only when that qubit
is in state |1⟩. The default `control=-2` adds no control.
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
    hilbert(N::Int, mat::AbstractMatrix, qubit::Int; control=-2)

Embed a 2×2 matrix on `qubit` in an `N`-qubit register, returning a sparse
matrix. Qubit 1 is the most significant tensor factor. Set `control` to a
distinct qubit to apply `mat` only when the control is |1⟩; the default adds
no control.
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

"""
    relabel_swap(ops::Vector{<:QuantumOps}, qubit_mapping::Vector{Int}=Int[])

Remove unconditional SWAP gates by relabeling subsequent operations. Return
`(relabeled_ops, final_mapping)`, where `final_mapping[i]` gives the final
physical label of qubit `i` in the relabeled circuit. Account for this mapping
when interpreting the final state or measurements.

The initial mapping is the identity unless supplied as a permutation of `1:N`;
newly encountered higher-index qubits extend it with identity entries. Controlled
SWAPs cannot be removed this way. An opaque whole-state `OpF` cannot be relabeled
after a permutation; use a block containing explicit operations instead.

```julia
ops = [Op("SWAP", 1, 2), Op("X", 1)]
relabeled_ops, final_mapping = relabel_swap(ops)
@assert length(relabeled_ops) == 1
@assert only(relabeled_ops).qubit == 2
@assert final_mapping == [2, 1]
```
"""
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
    apply_noise(state, op::QuantumOps, noise::NoiseModel)

Apply the model's one- or two-qubit channel to the operands of `op`, provided
the operation has `noisy=true`. This function applies only the noise channel;
`apply(state, op; noise)` applies the gate followed by its noise.
Statevectors and MPS sample a Kraus branch; density matrices sum all branches.
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
    apply(ops; N=get_N(ops), backend=:statevector, noise=false,
          track_measurements=false, layout=nothing, kwargs...)

Apply operations in their listed order to an existing state (`:statevector`,
`:density_matrix`, `:tensor`/`:mps`, or `CAMPSState`), or starting from `|0…0⟩`
in the requested `backend` (`:statevector`, `:density_matrix`, `:tensor`,
`:camps`, or `:phasepoly`).

By default, density matrices sum measurement and noise branches, while
statevectors, MPS, and CAMPS follow one sampled trajectory. Pass
`track_measurements=true` to return `(final_state, bits)`, conditioning density
matrices on the sampled measurement outcomes and recording each explicit
mid-circuit measurement bit (`MX`, `MY`, `MZ`, `MR`, `ifOp`) in operation
order. No routing or scheduling is performed unless `layout` is supplied.
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
