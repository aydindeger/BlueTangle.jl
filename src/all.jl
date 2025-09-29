import SparseArrays as sa
import LinearAlgebra as la
import ITensors as its
import ITensorMPS as it
import StatsBase as sb
import Pkg

using ForwardDiff, Optimisers, OptimKit

const AbstractVectorS = Union{AbstractVector, sa.SparseVector}
const AbstractMatrixS = Union{AbstractMatrix, sa.SparseMatrixCSC}

include("struct.jl")
include("decompose.jl")
include("hilbert.jl")
include("func.jl")
include("ops.jl")
include("qec.jl")
include("bit.jl")
include("layout.jl")
include("gates.jl")
include("noise.jl")
include("tensor.jl")
include("trotter.jl")
include("vqa.jl")
include("linalg.jl")
include("qasm.jl")

import Base: *
*(o::QuantumOps, state::AbstractVectorS) = apply(state,o)
*(o::QuantumOps, psi::it.MPS) = apply(psi,o)
*(tensor::its.ITensor,psi::it.MPS) = it.apply(tensor,psi) #exact

*(str::Tuple, state::AbstractVectorS) = apply(state,Op(str...))
*(vector_str::AbstractVector{<:Tuple}, state::AbstractVectorS) = apply(vector_str, state)

*(vector_str::Tuple, psi::it.MPS) =  apply(vector_str, psi)
*(vector_str::AbstractVector{<:Tuple}, psi::it.MPS) =  apply(vector_str, psi)

"""
    ⊗(x::AbstractMatrixS,y::AbstractMatrixS)=kron
"""
⊗(x::AbstractMatrixS,y::AbstractMatrixS)=foldl(kron,[x,y])

fields(m)=fieldnames(typeof(m))
attributes(m)=fields(m)

"""
ro3(x)=round(x,sigdigits=3)
"""
ro3(x)=round(x,sigdigits=3)

"""
ro10(x)=round(x,digits=10)
"""
ro10(x)=round(x,digits=10)