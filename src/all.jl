import SparseArrays as sa
import LinearAlgebra as la
import ITensors as it
import StatsBase as sb
import Pkg

using ForwardDiff, Optimisers, OptimKit

include("struct.jl")
include("decompose.jl")
include("hilbert.jl")
include("func.jl")
include("ops.jl")
include("bit.jl")
include("layout.jl")
include("gates.jl")
include("noise.jl")
include("tensor.jl")
include("trotter.jl")
include("vqe.jl")
include("linalg.jl")

import Base: *
*(o::QuantumOps, state::sa.SparseVector) = apply(state,o)
*(o::QuantumOps, psi::it.MPS) = apply(psi,o)
*(tensor::it.ITensor,psi::it.MPS) = it.apply(tensor,psi) #exact

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