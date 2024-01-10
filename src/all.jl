using SparseArrays, LinearAlgebra, Pkg
import StatsBase

include("struct.jl")
include("func.jl")
include("bit.jl")
include("gates.jl")
include("noise.jl")

using PyPlot

include("plot.jl")
include("trotter.jl")