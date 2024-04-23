
##
##========== expectation and correlations from measurement ==========

"""
`correlation(m::Measurement, qubits::Vector) -> Any`

Calculate correlation from a measurement on specified qubits.

- `m`: A `Measurement` object.
- `qubits`: Vector of qubits for measuring correlations, e.g. eg: qubits=[1,3] measure <Z1Z3>

Returns the calculated correlation.
"""
correlation(m::Measurement,qubits::Vector)=_sample_to_expectation((fock_basis.(m.bitstr,m.number_of_qubits),m.sample),qubits)

"""
Alias:
```
expect(m::Measurement, qubit::Int) -> Float64
expect(m::Measurement) -> Vector
```

Calculate the expectation value from a measurement.

- `m`: A `Measurement` object.
- `qubit`: (Optional) Specific qubit for which to calculate the expectation value.

Returns the expectation value(s).
"""
expect(m::Measurement,qubit::Int)=_sample_to_expectation((fock_basis.(m.bitstr,m.number_of_qubits),m.sample),[qubit])

##========== expectation and correlations from measurement ==========

##========== expectation from state & rho and operator ==========

"""
Alias for density matrix:
```
expect(state::sa.SparseVector, op::QuantumOps) -> Float64
expect(state::sa.SparseVector, op_str::String, qubit::Int) -> Float64
expect(state::sa.SparseVector, matrix::sa.SparseMatrixCSC) -> Float64
expect(state::sa.SparseVector, op_str::String) -> Vector{Float64}
```

Alias for state vector:
```
expect(rho::sa.SparseMatrixCSC, op::QuantumOps) -> Float64
expect(rho::sa.SparseMatrixCSC, op_str::String, qubit::Int) -> Float64
expect(rho::sa.SparseMatrixCSC, matrix::sa.SparseMatrixCSC) -> Float64
expect(rho::sa.SparseMatrixCSC, op_str::String) -> Vector{Float64}
```

Calculate the expectation value for quantum states or density matrices given an operator. This function has several forms depending on the input parameters:

- `expect(state::sa.SparseVector, op::QuantumOps)`: Computes the expectation value for a quantum state vector with a given operator.

- `expect(rho::sa.SparseMatrixCSC, op::QuantumOps)`: Computes the expectation value for a density matrix with a given operator.

- `expect(state::sa.SparseVector, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a quantum state vector with an operator specified as a string.

- `expect(rho::sa.SparseMatrixCSC, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a density matrix with an operator specified as a string.

- `expect(state::sa.SparseVector, op_str::String)`: Computes the expectation values for all qubits in a quantum state vector given an operator as a string.

- `expect(rho::sa.SparseMatrixCSC, op_str::String)`: Computes the expectation values for all qubits in a density matrix given an operator as a string.

- `expect(state::sa.SparseVector, matrix::sa.SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a state vector.

- `expect(rho::sa.SparseMatrixCSC, matrix::sa.SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a density matrix.

# Arguments
- `state::sa.SparseVector`: The quantum state vector.
- `rho::sa.SparseMatrixCSC`: The density matrix.
- `op::QuantumOps`: The quantum operator.
- `op_str::String`: The string representation of the operator.
- `qubit::Int`: The specific qubit index.
- `matrix::sa.SparseMatrixCSC`: The sparse matrix representation of the operator.

# Returns
- The expectation value as a `Float64` or a vector of `Float64` for multiple qubits.
"""
expect(state::sa.SparseVector,op::QuantumOps)=real(state'*op.expand(get_N(state))*state)
expect(rho::sa.SparseMatrixCSC,op::QuantumOps)=real(la.tr(rho*op.expand(get_N(rho))))

# expect(state::sa.SparseVector,op_str::String,qubit::Int)=real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state)
# expect(rho::sa.SparseMatrixCSC,op_str::String,qubit::Int)=real(la.tr(rho*expand_multi_op(op_str,[qubit],get_N(rho))))

expect(state::sa.SparseVector,op_str::String)=[real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state) for qubit=1:get_N(state)]
expect(rho::sa.SparseMatrixCSC,op_str::String)=[real(la.tr(rho*expand_multi_op(op_str,[qubit],get_N(rho)))) for qubit=1:get_N(rho)]

expect(state::sa.SparseVector,matrix::sa.SparseMatrixCSC)=real(state'*matrix*state)
expect(rho::sa.SparseMatrixCSC,matrix::sa.SparseMatrixCSC)=real(la.tr(rho*matrix))

"""
Alias:
```
correlation(state::sa.SparseVector, list_of_operators::String, qubits_applied::Vector) -> Float64
correlation(rho::sa.SparseMatrixCSC, list_of_operators::String, qubits_applied::Vector) -> Float64
```

Calculate the correlation for a given set of operators applied to specific qubits in either a quantum state vector or a density matrix. This function has two primary forms:

- `correlation(state::sa.SparseVector, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a quantum state vector (`state`) with a specified list of operators and qubits.

- `correlation(rho::sa.SparseMatrixCSC, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a density matrix (`rho`) with a specified list of operators and qubits.

The `corr_from_rho` function is an alias to `get_corr` for density matrices.

# Arguments
- `state::sa.SparseVector`: The quantum state vector.
- `rho::sa.SparseMatrixCSC`: The density matrix.
- `list_of_operators::String`: A string representing a list of operators, e.g., "Z,Z".
- `qubits_applied::Vector`: A vector of qubit indices on which the operators are applied.

# Returns
- `Float64`: The computed correlation value.

# Examples
```julia
# For a state vector
state = sa.SparseVector([...]) # define your state vector
correlation = correlation(state, "Z,Z", [2, 4])

# For a density matrix
rho = sa.SparseMatrixCSC([...]) # define your density matrix
correlation = correlation(rho, "Z,Z", [2, 4])
```

"""
function correlation(state::sa.SparseVector,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(state))
    return real(state'*matrix*state)
end

function correlation(rho::sa.SparseMatrixCSC,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(rho))
    return real(la.tr(rho*matrix))
end

##========== expectation from state & rho and operator ==========

#todo: is this expectation or correlation?
"""
`_sample_to_expectation(bprob::Tuple, qubits::Vector) -> Float64`

Calculates the expectation value for a specified set of qubits based on basis probability data from quantum measurements.

This function is designed to work with realistic quantum measurement data, where `bprob` is a tuple containing the basis states and their corresponding probabilities. The basis states are assumed to be in the Fock basis.

- `bprob`: A tuple where the first element is an array of basis states (in Fock basis) and the second element is an array of corresponding probabilities.
- `qubits`: A vector of qubit indices for which the expectation value is to be calculated.

The expectation value is computed by summing the probabilities of basis states where the specified qubits are in even parity and subtracting the probabilities of states with odd parity.

Returns the computed expectation value as a float.
"""
function _sample_to_expectation(bprob::Tuple,qubits::Vector)

focks=bprob[1]
prob=bprob[2]

probsum=0
for (f,p)=zip(focks,prob)
    if iseven(sum(f[qubits]))
        probsum+=p
    else
        probsum-=p
    end
end

return probsum

end

##========== expectation from measurement ==========

mag_moments_from_rho(rho::sa.SparseMatrixCSC,op_str::String;max_moment=12)=[mag_moments_from_rho(rho,op_str,ord) for ord=1:max_moment]

"""
`mag_moments_from_rho(rho::sa.SparseMatrixCSC, op_str::String, moment_order::Int) -> Float64`

Calculates the magnetization moments from a given density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `moment_order`: The order of the magnetization moment to compute.

Returns the computed magnetization moment of the specified order.
"""
function mag_moments_from_rho(rho::sa.SparseMatrixCSC,op_str::String,moment_order::Int)
    throw("fix")
    N=get_N(rho)
    mag_op=sum(hilbert_op(Op(op_str,i),N) for i=1:N)
    return real(la.tr(mag_op^moment_order*rho))
end

"""
`mag_cumulants_from_rho(rho::sa.SparseMatrixCSC, op_str::String, cumulant_order::Int; max_moment::Int=12) -> Float64`

Calculates the magnetization cumulants from a given density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `cumulant_order`: The order of the cumulant to compute.
- `max_moment`: (Optional) Maximum order of moments to consider.

Returns the computed magnetization cumulant of the specified order.
"""
function mag_cumulants_from_rho(rho::sa.SparseMatrixCSC,op_str::String,cumulant_order::Int;max_moment=12)
    moments=[mag_moments_from_rho(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments,cumulant_order)
end

function mag_cumulants_from_rho(rho::sa.SparseMatrixCSC,op_str::String;max_moment=12)
    moments=[mag_moments_from_rho(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments)
end

"""
calculates average magnetization moments from sample
"""
get_mag_moments(m::Measurement,moment_order::Int)=get_mag_moments(m.number_of_qubits,m.bitstr,m.sample,moment_order)
get_mag_moments(m::Measurement;max_moment=12)=[get_mag_moments(m::Measurement,ord) for ord=1:max_moment]

function get_mag_moments(number_of_qubits::Int,bitstr::Union{UnitRange,Vector},sample::Union{sa.SparseVector,Vector},moment_order::Int)
    mag=sum.(BlueTangle.mag_basis.(bitstr,number_of_qubits))
    return sum(mag.^(moment_order) .* sample)
end

##

cumulants_from_moments(moments::Vector)=[cumulants_from_moments(moments,n) for n=1:length(moments)]

function cumulants_from_moments(moments::Vector,n::Int)
    if n == 1
        return moments[1]
    else
        cumulant_n = moments[n]
        for m in 1:(n-1)
            cumulant_n -= binomial(n-1, m-1) * cumulants_from_moments(moments,m) * moments[n-m]
        end
        return cumulant_n
    end
end

function mag_cumulants_from_measurement(m::Measurement,cumulant_order::Int;max_moment=12)
    moments=[get_mag_moments(m,i) for i=1:max_moment]
    return cumulants_from_moments(moments,cumulant_order)
end


"""
MPS high order cumulants
"""
function calculate_high_moments(psi::it.MPS,sites::Vector,direction::String="Z",max_moment_order=10)

    N=length(psi)
    os = it.OpSum()

    for j=1:N
        os += direction,j
    end
    
    W0 = it.MPO(os,sites)

    moments_MPO=[]

    W=deepcopy(W0)
    for i=1:max_moment_order
        push!(moments_MPO,inner(psi',W,psi))
        W=it.apply(W,W0)
    end
    return moments_MPO

end


##
##========== Entanglement entropy ==========

"""
    entanglement_entropy(psi::sa.SparseVector) -> Float64

Calculates the entanglement entropy of a quantum state.

- `psi`: The quantum state vector (sa.SparseVector) for which to calculate the entropy.

Returns the calculated entanglement entropy.
"""
function entanglement_entropy(psi::sa.SparseVector)

    N=get_N(psi)
    partA=N ÷ 2
    
    mat=reshape(psi,(2^(partA),2^(N-partA)));
    
    spec=la.svdvals(Matrix(mat)).^2
    
    spec=spec[spec .> 0]
    
    return sum(-spec.*log.(spec)),-log.(spec)
    
end

"""
 partial_trace(rho) -> Matrix

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


function partial_trace(N::Int,state::sa.SparseVector,keep_index1::Int)

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

function partial_trace(N::Int,state::sa.SparseVector,keep_index1::Int,keep_index2::Int)

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


# N=6
# state=sa.normalize(sa.sparse(rand(ComplexF64,2^N)));


# keep_index1=4
# keep_index2=5
# traceout_ind=setdiff(1:N,[keep_index1,keep_index2])
# a=BlueTangle.partial_trace(state*state',fill(2,N),traceout_ind)

# b=partial_trace(N,state,keep_index1,keep_index2)

# println("two qubit implementation=",a ≈ b)


# traceout_ind=setdiff(1:N,[keep_index1])
# c=BlueTangle.partial_trace(state*state',fill(2,N),traceout_ind)
# d=partial_trace(N,state,keep_index1)

# println("one qubit implementation=",c ≈ d)


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
    entanglement_entropy(rho::sa.SparseMatrixCSC) -> Float64

Calculates the entanglement entropy of a density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) for which to calculate the entropy.

Returns the calculated entanglement entropy.
"""
function entanglement_entropy(rho::sa.SparseMatrixCSC)
    mat=bipartition_trace(rho)#trace half
    spec=la.svdvals(mat)
    spec=spec[spec .> 0]
    return sum(-spec.*log.(spec)),-log.(spec) 
end

### linear fit

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

ro3(x)=round(x,sigdigits=3)