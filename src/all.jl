import SparseArrays as sa
import LinearAlgebra as la
import ITensors as its
import ITensorMPS as it
import StatsBase as sb
import Serialization

using ForwardDiff, Optimisers, OptimKit
import PRIMA as pr

const AbstractVectorS = Union{AbstractVector,sa.SparseVector}
const AbstractMatrixS = Union{AbstractMatrix,sa.SparseMatrixCSC}

import Base: *, run

include("struct.jl")
include("decompose.jl")
include("hilbert.jl")
include("func.jl")
include("ops.jl")
include("qec.jl")
include("qec_steane.jl")
include("qec_round.jl")
include("qec_faults.jl")
include("qec_shor.jl")
include("qec_prepare.jl")
include("qec_magic.jl")
include("qec_codes.jl")
include("qec_distance.jl")
include("qec_distillation.jl")
include("qec_factory.jl")
include("qec_iswap.jl")
include("bit.jl")
include("layout.jl")
include("gates.jl")
include("noise.jl")
include("backends/statevector.jl")
include("backends/density_matrix.jl")
include("backends/tensor.jl")
include("trotter.jl")
include("variational_state.jl")
include("vqa.jl")
include("linalg.jl")
include("qasm.jl")
include("plotqt.jl")

"""
Displaying docstrings in VSCode
https://discourse.julialang.org/t/displaying-docstrings-in-notebook-in-vscode-errors/107627/4
"""
macro docVSC(fn)
    return :(display("text/markdown", Base.doc($fn)))
end

*(o::QuantumOps, state::AbstractVectorS) = apply(state, o)
*(o::QuantumOps, rho::AbstractMatrixS) = apply(rho, o)
*(o::QuantumOps, psi::it.MPS) = apply(psi, o)
*(tensor::its.ITensor, psi::it.MPS) = it.apply(tensor, psi) #exact

*(str::Tuple, state::AbstractVectorS) = apply(state, Op(str...))
*(vector_str::AbstractVector{<:Tuple}, state::AbstractVectorS) = apply(state, vector_str)

*(str::Tuple, psi::it.MPS) = apply(psi, Op(str...))
*(vector_str::AbstractVector{<:Tuple}, psi::it.MPS) = apply(psi, vector_str)

"""
    ⊗(x::AbstractMatrixS,y::AbstractMatrixS)=kron
"""
⊗(x::AbstractMatrixS, y::AbstractMatrixS) = kron(x, y)

fields(m) = fieldnames(typeof(m))

"""
ro3(x)=round(x,sigdigits=3)
"""
ro3(x) = round(x, sigdigits=3)

"""
ro10(x)=round(x,digits=10)
"""
ro10(x) = round(x, digits=10)
