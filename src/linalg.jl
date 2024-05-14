function la.ishermitian(op::Op)
    if op.mat isa AbstractMatrix && la.ishermitian(op.mat)
        return true
    else
        return false
    end
end

Base.adjoint(op::OpF)=op
Base.adjoint(op::ifOp)=op

function Base.adjoint(op::Op)

    matf = op.mat
    name = op.name * "†"

    if la.ishermitian(op)
        return op #return original op
    end

    if matf isa Function
        adjmatf = x -> adjoint(matf(x))
    elseif matf isa la.Adjoint
        name = op.name[1:end-1]
        adjmatf = matf.parent
    else
        adjmatf = adjoint(matf)
    end
    return Op(name, adjmatf, op.qubit, op.target_qubit; type=op.type, noisy=op.noisy, control=op.control)

end

function Base.adjoint(ops::Vector{T}) where T<:QuantumOps
    return reverse(Base.adjoint.(ops))
end

function Base.adjoint(circ::Circuit)
    circ2=deepcopy(circ)
    for ops=circ2.layers
        replace!(op -> la.ishermitian(op) ? op : adjoint(op), ops)
    end

    reverse!(circ2.layers)

    return circ2
end

"""
    isunitary(mat::AbstractMatrix) --> Bool
"""
isunitary(mat::AbstractMatrix)=isapprox(mat'mat,Matrix(la.I, size(mat,1),size(mat,1)),atol=1e-12)

"""
    ishermitian(mat::AbstractMatrix) --> Bool
"""
ishermitian(mat::AbstractMatrix)=la.ishermitian(mat)

"""
    sparsevector(vec::AbstractArray)=sa.SparseVector(vec)
"""
sparsevector(vec::AbstractArray)=sa.SparseVector(vec)