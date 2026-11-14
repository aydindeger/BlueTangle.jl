function la.ishermitian(op::Op)
    if op.mat isa AbstractMatrix && la.ishermitian(op.mat)
        return true
    else
        return false
    end
end

Base.adjoint(op::OpF)=op
Base.adjoint(op::ifOp)=op

# Preserve the original function's parameter count when taking a variational
# gate's adjoint, so VQA can allocate the same parameter slice.
struct _AdjointGate{F} <: Function
    parent::F
end
(f::_AdjointGate)(args...) = adjoint(f.parent(args...))
_find_argument_number(f::_AdjointGate) = _find_argument_number(f.parent)

function Base.adjoint(op::Op)

    matf = op.mat
    name = op.name * "†"

    if la.ishermitian(op)
        return op #return original op
    end

    if matf isa _AdjointGate
        name = op.name[1:end-1]
        adjmatf = matf.parent
    elseif matf isa Function
        adjmatf = _AdjointGate(matf)
    elseif matf isa la.Adjoint
        name = op.name[1:end-1]
        adjmatf = matf.parent
    else
        adjmatf = adjoint(matf)
    end
    return Op(name, adjmatf, op.qubit, op.target_qubit; noisy=op.noisy, control=op.control)

end

function Base.adjoint(ops::Vector{T}) where T<:QuantumOps
    return reverse(Base.adjoint.(ops))
end

"""
    isunitary(mat::AbstractMatrix) --> Bool
"""
isunitary(mat::AbstractMatrix)=isapprox(mat'mat,Matrix(la.I, size(mat,1),size(mat,1)),atol=1e-12)

"""
    sparsevector(vec::AbstractArray)=sa.SparseVector(vec)
"""
sparsevector(vec::AbstractArray)=sa.SparseVector(vec)

##


"""
    same_up_to_global_phase(psi::AbstractVectorS, phi::AbstractVectorS)

checks if following true:
psi==exp(im*θ)*phi
"""
function same_up_to_global_phase(psi::AbstractVectorS, phi::AbstractVectorS)
    inner_product = la.dot(phi,psi)
    global_phase = la.angle(inner_product)
    is_same = isapprox(psi, exp(im * global_phase) * phi)
    return is_same, global_phase
end


## 
##========== partial trace  ==========

# Basis offsets for selected subsystems, with subsystem 1 the most significant.
# Compute these once instead of rebuilding full indices for every matrix element.
function _partial_trace_offsets(dims, subsystems)
    sub_dims = dims[subsystems]
    strides = [prod(dims[i+1:end]) for i in subsystems]
    offsets = zeros(Int, prod(sub_dims))
    for i in eachindex(offsets)
        value = i - 1
        for k in reverse(eachindex(sub_dims))
            value, digit = divrem(value, sub_dims[k])
            offsets[i] += digit * strides[k]
        end
    end
    return offsets
end

"""
    partial_trace(state::AbstractVectorS, keep_qubits::AbstractVector)

Return the reduced density matrix of the selected qubits, in ascending qubit
order. Nonlocal qubits are supported. Contract the amplitudes directly, without
constructing the full density matrix of the state.
"""
function partial_trace(state::AbstractVectorS, keep_qubits::AbstractVector)
    N = get_N(state)
    dims = fill(2, N)
    trace_out = setdiff(1:N, keep_qubits)
    keep = setdiff(1:N, trace_out)
    keep_offsets = _partial_trace_offsets(dims, keep)
    trace_offsets = _partial_trace_offsets(dims, trace_out)
    reduced_dim = length(keep_offsets)
    ρ_reduced = zeros(ComplexF64, reduced_dim, reduced_dim)

    for j in 1:reduced_dim, i in 1:reduced_dim
        value = 0.0 + 0.0im
        for offset in trace_offsets
            value += state[keep_offsets[i] + offset + 1] *
                     conj(state[keep_offsets[j] + offset + 1])
        end
        ρ_reduced[i, j] = value
    end
    return ρ_reduced
end

"""
    partial_trace(ρ::AbstractMatrixS, dims::Vector, trace_out::Vector)

Trace out the specified subsystems of a density matrix whose subsystem
sizes are `dims`. The remaining subsystems retain their original order.
"""
function partial_trace(ρ::AbstractMatrixS, dims::Vector, trace_out::Vector)
    all(d -> d isa Integer && d > 0, dims) ||
        throw(ArgumentError("subsystem dimensions must be positive integers"))
    all(q -> q isa Integer && 1 <= q <= length(dims), trace_out) && allunique(trace_out) ||
        throw(ArgumentError("trace_out must contain distinct valid subsystem indices"))
    total_dim = prod(dims)
    size(ρ) == (total_dim, total_dim) ||
        throw(DimensionMismatch("density matrix size must match the product of dims"))
    keep = setdiff(1:length(dims), trace_out)
    keep_offsets = _partial_trace_offsets(dims, keep)
    trace_offsets = _partial_trace_offsets(dims, sort(trace_out))
    reduced_dim = length(keep_offsets)
    ρ_reduced = zeros(ComplexF64, reduced_dim, reduced_dim)

    for j in 1:reduced_dim, i in 1:reduced_dim
        value = 0.0 + 0.0im
        for offset in trace_offsets
            value += ρ[keep_offsets[i] + offset + 1, keep_offsets[j] + offset + 1]
        end
        ρ_reduced[i, j] = value
    end
    return ρ_reduced
end

"""
    bipartition_trace(rho) -> Matrix

Computes the partial trace of a density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) for which to perform the partial trace.

Returns the resulting matrix after performing the partial trace.
"""
function bipartition_trace(rho::AbstractMatrixS)
    dim=Int(get_N(rho)/2)
    d=2^dim
    ptr = zeros(Complex{Float64}, d, d)
    for i = 1:d, j = 1:d
        for k = 0:d-1
            ptr[i, j] += rho[i+k*d, j+k*d]
        end
    end
    return ptr
end


"""
    partial_trace(state::AbstractVectorS,keep_index1::Int)
"""
function partial_trace(state::AbstractVectorS,keep_index1::Int)

    N=get_N(state)

    qubit_index =  N - keep_index1 + 1

    size_left = 2^(qubit_index - 1)
    size_middle = 2
    size_right = 2^(N - qubit_index)

    # Reshape the state vector to correctly isolate the two-qubit subsystem
    tensor = reshape(state, (size_left, size_middle, size_right))

    # Initialize the reduced density matrix for the qubit of interest
    ρ_reduced = zeros(ComplexF64, size_middle, size_middle)

    # Sum over the other qubits to obtain the reduced density matrix
    for i in 1:size_left
        for j in 1:size_right
            ρ_reduced += tensor[i, :, j] * tensor[i, :, j]'
        end
    end

    ρ_reduced

end

"""
    partial_trace(state::AbstractVectorS,keep_index1::Int,keep_index2::Int)
    index1 and index2 are local
"""
function partial_trace(state::AbstractVectorS,keep_index1::Int,keep_index2::Int)

    N=get_N(state)

    if abs(keep_index1-keep_index2)>1
        throw("must be local")
    end

    qubit_index1 = N-minimum([keep_index1,keep_index2])#this fixes enumeration.
    qubit_index2 = qubit_index1+1

    size_left = 2^(qubit_index1 - 1)
    size_middle = 4  # Since we're dealing with a two-qubit subsystem, this should directly be 2^2
    size_right = 2^(N - qubit_index2)

    # Reshape the state vector to correctly isolate the two-qubit subsystem
    tensor = reshape(state, (size_left, size_middle, size_right))

    # Initialize the reduced density matrix for the two-qubit subsystem as a 4x4 matrix
    ρ_reduced = zeros(ComplexF64, 4, 4)

    # Adjust the loop to correctly compute the reduced density matrix
    for i in 1:size_left
        for k in 1:size_right
            # Isolate the two-qubit subsystem
            mat = tensor[i, :, k]  # This should be a 4-element vector for the two-qubit subsystem
            ρ_reduced += mat * mat'
        end
    end

    return ρ_reduced

end


# Define a function to perform the partial trace over the second subsystem
function partial_trace_second_subsystem(ρ, dimA, dimB)
    # ρ is the full density matrix of the system A ⊗ B
    # dimA is the dimension of subsystem A
    # dimB is the dimension of subsystem B
    
    # Initialize the reduced density matrix for subsystem A
    ρA = zeros(ComplexF64, dimA, dimA)
    
    # Perform the partial trace over subsystem B
    for i in 1:dimA
        for j in 1:dimA
            for k in 1:dimB
                ρA[i, j] += ρ[(i-1)*dimB + k, (j-1)*dimB + k]
            end
        end
    end
    
    return ρA
end

"""
    stinespring_dilation(kraus_ops::Vector{AbstractMatrix}) -> Unitary Matrix
"""
function stinespring_dilation(kraus_ops::Vector{<:AbstractMatrix})
    n = size(kraus_ops[1], 1)  # Dimension of the system
    d = n * length(kraus_ops) 
    Q = vcat(kraus_ops...)

    if d > n
        U = hcat(Q, la.nullspace(Q')) # Complete the orthonormal basis to form a unitary matrix
    else
        U = Q
    end
    return isunitary(U) ? U : throw("not unitary. something's wrong")
end

##========== partial trace  ==========

"""
    reduced_row_echelon(generator::AbstractMatrix)
"""
function reduced_row_echelon(generator::AbstractMatrix{Int})
    num_rows, num_cols = size(generator)
    transform_rows = Matrix{Int}(la.I(num_rows))
    pivots = []
    i = 1  # row index
    j = 1  # column index

    while i <= num_rows && j <= num_cols
        # Find the pivot in column j starting from row i
        pivot_row = findfirst(generator[i:end, j] .!= 0)
        if pivot_row !== nothing
            pivot_row += i - 1  # Adjust index because of slicing
            # Swap rows if necessary
            if pivot_row != i
                generator[[i, pivot_row], :] = generator[[pivot_row, i], :]
                transform_rows[[i, pivot_row], :] = transform_rows[[pivot_row, i], :]
            end

            # Record pivot position
            push!(pivots, (i, j))

            # Eliminate entries in other rows
            for k in 1:num_rows
                if k != i && generator[k, j] == 1
                    generator[k, :] .= mod.(generator[k, :] .+ generator[i, :], 2)
                    transform_rows[k, :] .= mod.(transform_rows[k, :] .+ transform_rows[i, :], 2)
                end
            end
            i += 1
        end
        j += 1  # Always move to the next column
    end

    rank = length(pivots)

    # Reorder columns to place pivot columns first
    pivot_cols = [j for (i, j) in pivots]
    non_pivot_cols = setdiff(1:num_cols, pivot_cols)
    permuted_indices = vcat(pivot_cols, non_pivot_cols)

    # Reorder generator and transform_cols accordingly
    generator = generator[:, permuted_indices]
    transform_cols = Matrix{Int}(la.I(num_cols))[:, permuted_indices]

    return generator .% 2, rank, transform_rows .% 2, transform_cols
end

function rank_of_rref(A::AbstractMatrix{Int})
    len_cols=size(A,2)÷2
    non_zero_rows = sum(any(A[i, 1:len_cols] .== 1) for i in 1:size(A, 1))
    return non_zero_rows
end

##
##========== fit  ==========


"""
Simple linear model: y = a + b*x
"""
_linear_model(x, a, b)=a .+ b .* x
_quadratic_model(x,a,b,c)=a .* x .^2 .+ b .*x .+ c

"""
    linear_fit(xdata::Vector, ydata::Vector)

Compute the coefficients a and b for the linear fit of the given data.

# Arguments
- `xdata`: Array of x values.
- `ydata`: Array of corresponding y values.

# Returns
- `a`: Intercept of the linear fit.
- `b`: Slope of the linear fit.
- `se_a`: Standard error of the intercept.
- `se_b`: Standard error of the slope.

# Description
This function fits a simple linear model (y = a + b*x) to the provided data points.
It checks if the lengths of xdata and ydata are the same and then calculates the coefficients
for the linear fit. Additionally, it computes the standard errors for both the slope and the intercept.

The function uses the least squares method for the linear regression. The standard errors are calculated
based on the residual sum of squares and the total sum of squares for the x values.
"""
function linear_fit(xdata::Vector, ydata::Vector)
    n = length(xdata)
    n == length(ydata) || throw(DimensionMismatch("xdata and ydata must have the same length"))
    n >= 2 || throw(ArgumentError("a line requires at least two points"))

    # Center before summing products to avoid subtracting nearly equal large
    # numbers when the x coordinates have a large common offset.
    x_mean = sum(xdata) / n
    y_mean = sum(ydata) / n
    sst_x = sum((x - x_mean)^2 for x in xdata)
    sst_x > 0 || throw(ArgumentError("xdata must contain at least two distinct values"))
    b = sum((x - x_mean) * (y - y_mean) for (x, y) in zip(xdata, ydata)) / sst_x
    a = y_mean - b * x_mean
    rss = sum((y - (a + b * x))^2 for (x, y) in zip(xdata, ydata))
    # Two points fix a line but cannot estimate residual uncertainty.
    residual_variance = n == 2 ? NaN : rss / (n - 2)
    se_b = sqrt(residual_variance / sst_x)
    se_a = sqrt(residual_variance * (1 / n + x_mean^2 / sst_x))
    return a, b, se_a, se_b
end


"""
    quadratic_fit(xdata::Vector{Float64}, ydata::Vector{Float64}) -> Vector{Float64}

Fit a quadratic function to the given data.

# Arguments
- `xdata::Vector{Float64}`: A vector of x-coordinates.
- `ydata::Vector{Float64}`: A vector of y-coordinates corresponding to `xdata`.

# Returns
- `Vector{Float64}`: The coefficients `[a, b, c]` of the fitted quadratic function `y = ax^2 + bx + c`.

# Description
This function performs a least squares quadratic fit to the input data.
"""
function quadratic_fit(xdata::Vector{Float64}, ydata::Vector{Float64})
    if length(xdata) != length(ydata)
        error("Vectors xdata and ydata must have the same length")
    end

    n = length(xdata)
    X = ones(n, 3)
    for i in 1:n
        X[i, 1] = xdata[i]^2
        X[i, 2] = xdata[i]
    end
    Y = reshape(ydata, n, 1)  # Convert y to a column vector

    # QR least squares avoids squaring the condition number via X'X.
    beta = X \ Y

    return beta # the coefficients [a, b, c]
end


"""
 error_mitigate_data(xdata::Vector, ydata::Vector)

Perform error mitigation on a dataset by fitting a linear model and extracting the estimate and standard error.

This function takes two vectors `xdata` and `ydata` which represent the independent and dependent variables of a dataset, respectively.

# Arguments
- `xdata::Vector`: The independent variable data points.
- `ydata::Vector`: The dependent variable data points, corresponding to each xdata point.

# Returns
- `est`: The estimated intercept from the linear fit.
- `se`: The standard error of the estimated intercept.
- `fit_plot`: A tuple containing the x-values from 0 to the last element of `xdata` and the corresponding fitted y-values from the model.
"""
function error_mitigate_data(xdata::Vector,ydata::Vector)
    a, b, se_a, se_b = linear_fit(xdata, ydata)
    est=round(a,sigdigits=5)

    dx=(xdata[2]-xdata[1])/4
    fit_plot=(0:dx:xdata[end],[_linear_model(x,a,b) for x=0:dx:xdata[end]])
    est,se_a,fit_plot
end

##========== fit  ==========