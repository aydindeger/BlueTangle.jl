import SparseArrays as sa
import LinearAlgebra as la
import ITensors as it
import StatsBase as sb
import Pkg

using ForwardDiff, Optimisers, OptimKit

include("struct.jl")
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

import Base: *
*(o::QuantumOps, state::sa.SparseVector) = apply(o,state)