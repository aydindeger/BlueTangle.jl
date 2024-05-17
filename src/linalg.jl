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
    sparsevector(vec::AbstractArray)=sa.SparseVector(vec)
"""
sparsevector(vec::AbstractArray)=sa.SparseVector(vec)

## 
##========== partial trace  ==========

"""
    partial_trace(ρ::sa.SparseVector, dims::Vector, keep_qubits::Vector)

    This function can handle nonlocal operations.
"""
function partial_trace(state::sa.SparseVector, keep_qubits::Vector)
    N=get_N(state)
    return partial_trace(state*state', fill(2,N),setdiff(1:N,keep_qubits))
end

function partial_trace(ρ::sa.SparseMatrixCSC, dims::Vector, trace_out::Vector)
    # julia version of algorithm: https://copyprogramming.com/howto/how-to-take-partial-trace

    # Calculate the total dimension of the system and initialize the reduced density matrix
    total_dim = prod(dims)
    trace_dims = [dims[i] for i in trace_out]
    keep_dims = deleteat!(copy(dims), trace_out)
    reduced_dim = prod(keep_dims)
    ρ_reduced = zeros(ComplexF64, reduced_dim, reduced_dim)

    # Helper function to calculate the index in the reduced density matrix
    function calc_index(indices, dims)
        idx = 0
        for i in 1:length(indices)
            idx = idx * dims[i] + indices[i] - 1
        end
        return idx + 1
    end

    # Iterate over the elements of the reduced density matrix
    for i in 1:reduced_dim
        for j in 1:reduced_dim
            # Convert linear indices to multi-indices for the kept subsystems
            i_indices = [div((i-1), prod(keep_dims[k+1:end])) % keep_dims[k] + 1 for k in 1:length(keep_dims)]
            j_indices = [div((j-1), prod(keep_dims[k+1:end])) % keep_dims[k] + 1 for k in 1:length(keep_dims)]
            
            # Sum over the traced out dimensions
            sum = 0.0
            for k in 1:prod(trace_dims)
                # Convert linear index to multi-indices for the traced out subsystems
                k_indices = [div((k-1), prod(trace_dims[l+1:end])) % trace_dims[l] + 1 for l in 1:length(trace_dims)]
                
                # Combine indices and calculate linear indices for the full matrix
                full_i_indices = copy(i_indices)
                full_j_indices = copy(j_indices)
                insert!(full_i_indices, trace_out[1], k_indices[1])
                insert!(full_j_indices, trace_out[1], k_indices[1])
                for m in 2:length(trace_out)
                    insert!(full_i_indices, trace_out[m], k_indices[m])
                    insert!(full_j_indices, trace_out[m], k_indices[m])
                end
                
                lin_idx_i = calc_index(full_i_indices, dims)
                lin_idx_j = calc_index(full_j_indices, dims)
                
                sum += ρ[lin_idx_i, lin_idx_j]
            end
            ρ_reduced[i, j] = sum
        end
    end

    return ρ_reduced
end

"""
    bipartition_trace(rho) -> Matrix

Computes the partial trace of a density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) for which to perform the partial trace.

Returns the resulting matrix after performing the partial trace.
"""
function bipartition_trace(rho::sa.SparseMatrixCSC)
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


function partial_trace(state::sa.SparseVector,keep_index1::Int)

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
    partial_trace(state::sa.SparseVector,keep_index1::Int,keep_index2::Int)
    index1 and index2 are local
"""
function partial_trace(state::sa.SparseVector,keep_index1::Int,keep_index2::Int)

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

##========== partial trace  ==========

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

    if n != length(ydata)
        throw("xdata and ydata must have the same length")
    end

    # Calculate sums needed for the coefficients
    sum_x = sum(xdata)
    sum_y = sum(ydata)
    sum_x2 = sum(xdata .* xdata)
    sum_xy = sum(xdata .* ydata)

    # Calculate coefficients
    b = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)
    a = (sum_y - b * sum_x) / n

    rss = sum((y - (a + b * x))^2 for (x, y) in zip(xdata, ydata))

    # Calculate total sum of squares for x
    x_mean = sum_x / n
    sst_x = sum((x - x_mean)^2 for x in xdata)

    # Standard Error of the Slope (SE(b))
    se_b = sqrt(rss / (n - 2) / sst_x)

    # Standard Error of the Intercept (SE(a))
    se_a = se_b * sqrt(sum_x2 / n)

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

    beta = inv(X' * X) * X' * Y

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