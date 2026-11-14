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


function correlation(state::AbstractVectorS,qubits::Vector)
    
a,b=sample_exact(state)

return _sample_to_expectation((fock_basis.(a,get_N(state)),b),qubits)

end

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
expect(state::AbstractVectorS, op::QuantumOps) -> Float64
expect(state::AbstractVectorS, op_str::String, qubit::Int) -> Float64
expect(state::AbstractVectorS, matrix::sa.SparseMatrixCSC) -> Float64
expect(state::AbstractVectorS, op_str::String) -> Vector{Float64}
```

Alias for state vector:
```
expect(rho::AbstractMatrixS, op::QuantumOps) -> Float64
expect(rho::AbstractMatrixS, op_str::String, qubit::Int) -> Float64
expect(rho::AbstractMatrixS, matrix::sa.SparseMatrixCSC) -> Float64
expect(rho::AbstractMatrixS, op_str::String) -> Vector{Float64}
```

Calculate the expectation value for quantum states or density matrices given an operator. This function has several forms depending on the input parameters:

- `expect(state::AbstractVectorS, op::QuantumOps)`: Computes the expectation value for a quantum state vector with a given operator.

- `expect(rho::AbstractMatrixS, op::QuantumOps)`: Computes the expectation value for a density matrix with a given operator.

- `expect(state::AbstractVectorS, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a quantum state vector with an operator specified as a string.

- `expect(rho::AbstractMatrixS, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a density matrix with an operator specified as a string.

- `expect(state::AbstractVectorS, op_str::String)`: Computes the expectation values for all qubits in a quantum state vector given an operator as a string.

- `expect(rho::AbstractMatrixS, op_str::String)`: Computes the expectation values for all qubits in a density matrix given an operator as a string.

- `expect(state::AbstractVectorS, matrix::sa.SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a state vector.

- `expect(rho::AbstractMatrixS, matrix::sa.SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a density matrix.

# Arguments
- `state::AbstractVectorS`: The quantum state vector.
- `rho::AbstractMatrixS`: The density matrix.
- `op::QuantumOps`: The quantum operator.
- `op_str::String`: The string representation of the operator.
- `qubit::Int`: The specific qubit index.
- `matrix::sa.SparseMatrixCSC`: The sparse matrix representation of the operator.

# Returns
- The expectation value as a `Float64` or a vector of `Float64` for multiple qubits.
"""
expect(state::AbstractVectorS,op::QuantumOps)=real(state'*op.expand(get_N(state))*state)
expect(rho::AbstractMatrixS,op::QuantumOps)=real(la.tr(rho*op.expand(get_N(rho))))

# expect(state::AbstractVectorS,op_str::String,qubit::Int)=real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state)
# expect(rho::AbstractMatrixS,op_str::String,qubit::Int)=real(la.tr(rho*expand_multi_op(op_str,[qubit],get_N(rho))))

expect(state::AbstractVectorS,op_str::String)=[real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state) for qubit=1:get_N(state)]
expect(rho::AbstractMatrixS,op_str::String)=[real(la.tr(rho*expand_multi_op(op_str,[qubit],get_N(rho)))) for qubit=1:get_N(rho)]

expect(state::AbstractVectorS,matrix::AbstractMatrixS)=real(state'*matrix*state)
expect(rho::AbstractMatrixS,matrix::AbstractMatrixS)=real(la.tr(rho*matrix))

"""
Alias:
```
correlation(state::AbstractVectorS, list_of_operators::String, qubits_applied::Vector) -> Float64
correlation(rho::AbstractMatrixS, list_of_operators::String, qubits_applied::Vector) -> Float64
```

Calculate the correlation for a given set of operators applied to specific qubits in either a quantum state vector or a density matrix. This function has two primary forms:

- `correlation(state::AbstractVectorS, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a quantum state vector (`state`) with a specified list of operators and qubits.

- `correlation(rho::AbstractMatrixS, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a density matrix (`rho`) with a specified list of operators and qubits.

The `corr_from_rho` function is an alias to `get_corr` for density matrices.

# Arguments
- `state::AbstractVectorS`: The quantum state vector.
- `rho::AbstractMatrixS`: The density matrix.
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
function correlation(state::AbstractVectorS,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(state))
    return real(state'*matrix*state)
end

function correlation(rho::AbstractMatrixS,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(rho))
    return real(la.tr(rho*matrix))
end

##========== expectation from state & rho and operator ==========

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
    parity = false
    for q in qubits
        parity = xor(parity, isodd(f[q]))
    end
    if !parity
        probsum+=p
    else
        probsum-=p
    end
end

return probsum

end

##========== expectation from measurement ==========

mag_moments(rho::AbstractMatrixS,op_str::String;max_moment=12)=[mag_moments(rho,op_str,ord) for ord=1:max_moment]

"""
`mag_moments(rho::AbstractMatrixS, op_str::String, moment_order::Int) -> Float64`

Calculates the magnetization moments from a given density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `moment_order`: The order of the magnetization moment to compute.

Returns the computed magnetization moment of the specified order.
"""
function mag_moments(rho::AbstractMatrixS,op_str::String,moment_order::Int)
    moment_order >= 0 || throw(ArgumentError("moment order must be nonnegative"))
    N = get_N(rho)
    moment_order == 0 && return real(la.tr(rho))
    N == 0 && return zero(real(la.tr(rho)))
    mag_op = sum(Op(op_str, i).expand(N) for i in 1:N)
    return real(la.tr(mag_op^moment_order * rho))
end

"""
`mag_cumulants(rho::AbstractMatrixS, op_str::String, cumulant_order::Int; max_moment::Int=12) -> Float64`

Calculates the magnetization cumulants from a given density matrix.

- `rho`: The density matrix (sa.SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `cumulant_order`: The order of the cumulant to compute.
- `max_moment`: (Optional) Maximum order of moments to consider.

Returns the computed magnetization cumulant of the specified order.
"""
function mag_cumulants(rho::AbstractMatrixS,op_str::String,cumulant_order::Int;max_moment=12)
    moments=[mag_moments(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments,cumulant_order)
end

function mag_cumulants(rho::AbstractMatrixS,op_str::String;max_moment=12)
    moments=[mag_moments(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments)
end

"""
calculates average magnetization moments from sample
"""
mag_moments(m::Measurement,moment_order::Int)=mag_moments(m.number_of_qubits,m.bitstr,m.sample,moment_order)
mag_moments(m::Measurement;max_moment=12)=[mag_moments(m::Measurement,ord) for ord=1:max_moment]

function mag_moments(number_of_qubits::Int,bitstr::Union{UnitRange,Vector},sample::Union{sa.SparseVector,Vector},moment_order::Int)
    moment = 0.0
    for i in eachindex(bitstr, sample)
        # Read the same N low bits as mag_basis, without allocating bit vectors.
        ones = 0
        for shift in 0:number_of_qubits-1
            ones += (bitstr[i] >> shift) & 1
        end
        magnetization = number_of_qubits - 2ones
        # High moments can exceed Int even for modest qubit counts (40^12).
        moment += float(magnetization)^moment_order * sample[i]
    end
    return moment
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

function mag_cumulants(m::Measurement,cumulant_order::Int;max_moment=12)
    moments=[mag_moments(m,i) for i=1:max_moment]
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
        push!(moments_MPO, real(it.inner(psi', W, psi)))
        W=it.apply(W,W0)
    end
    return moments_MPO

end


##
##========== Entanglement entropy ==========

"""
    entanglement_entropy(psi::AbstractVectorS) -> Float64

Calculates the entanglement entropy of a quantum state.

- `psi`: The quantum state vector (sa.SparseVector) for which to calculate the entropy.

Returns `(entropy, spectrum)` by default for compatibility. Set
`spectrum_bool=false` to return only the entropy, or `true` for the tuple.
The spectrum contains `-log(p)` for nonzero Schmidt probabilities `p`.
"""
function entanglement_entropy(psi::AbstractVectorS; spectrum_bool::Union{Nothing,Bool}=nothing)
    N = get_N(psi)
    partA = N ÷ 2
    amplitudes = reshape(psi, (2^partA, 2^(N-partA)))
    probabilities = la.svdvals(Matrix(amplitudes)).^2
    probabilities = probabilities[probabilities .> 0]
    entropy = -sum(probabilities .* log.(probabilities))
    # Retain the historical tuple default; callers can explicitly request a scalar.
    return spectrum_bool === false ? entropy : (entropy, -log.(probabilities))
end

"""
    entanglement_entropy(rho::AbstractMatrix; spectrum_bool=nothing)

Return `(entropy, spectrum)` for the reduced density matrix of the last
`floor(N/2)` qubits, using natural logarithms. Set `spectrum_bool=false` for
just the entropy, or `true` for the tuple. The no-keyword tuple return is
retained for compatibility. For a mixed global state, reduced entropy includes
both quantum correlations and statistical mixing; it is not an entanglement
measure on its own.
"""
function entanglement_entropy(rho::AbstractMatrixS; spectrum_bool::Union{Nothing,Bool}=nothing)
    N = get_N(rho)
    reduced = partial_trace(rho, fill(2, N), collect(1:N-N÷2))
    probabilities = la.eigvals(la.Hermitian(Matrix(reduced)))
    minimum(probabilities) >= -1e-12 ||
        throw(ArgumentError("density matrix must be positive semidefinite"))
    probabilities = sort(filter(>(0), probabilities); rev=true)
    entropy = -sum(probabilities .* log.(probabilities))
    return spectrum_bool === false ? entropy : (entropy, -log.(probabilities))
end
