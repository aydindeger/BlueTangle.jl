##
##========== expectation and correlations from measurement ==========

"""
    correlation(m::Measurement, qubits::Vector)

Return the expectation of the product of measured eigenvalues on `qubits`.
Each recorded bit contributes `(-1)^bit`; for Z-basis data, `[1, 3]` therefore
gives ⟨Z₁Z₃⟩. The calculation uses the probabilities or signed weights in `m`.
"""
correlation(m::Measurement,qubits::Vector)=_sample_to_expectation((fock_basis.(m.bitstr,m.number_of_qubits),m.sample),qubits)


function correlation(state::AbstractVectorS,qubits::Vector)
    
a,b=sample_exact(state)

return _sample_to_expectation((fock_basis.(a,get_N(state)),b),qubits)

end

"""
    expect(m::Measurement, qubit::Int)

Return the expectation value on one qubit from the outcome weights. Use
`m.expect` for the vector of expectations on all qubits, in the measurement's
recorded basis.
"""
expect(m::Measurement,qubit::Int)=_sample_to_expectation((fock_basis.(m.bitstr,m.number_of_qubits),m.sample),[qubit])

##========== expectation and correlations from measurement ==========

##========== expectation from state & rho and operator ==========

"""
    expect(state, op::QuantumOps)
    expect(state, matrix::AbstractMatrix)
    expect(state, op_str::String)

Calculate expectations for a statevector or density matrix. An `Op` specifies
the qubits on which the observable acts; a matrix acts on the full register.
A single-qubit operator string returns one expectation per qubit. For one
selected qubit, use `expect(state, Op("Z", q))`. Results are real-valued.

```julia
state = apply(zero_state(2), Op("X", 2))
rho = to_rho(state)
@assert expect(state, "Z") ≈ [1, -1]
@assert expect(rho, Op("Z", 2)) ≈ -1
@assert expect(state, gate.Z ⊗ gate.Z) ≈ -1
measurement = measure(state)
@assert expect(measurement, 2) ≈ measurement.expect[2] ≈ -1
```
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
    correlation(state, list_of_operators::String, qubits_applied::Vector)

Return the expectation of a product of single-qubit observables in a statevector
or density matrix. Pair each comma-separated operator with the corresponding
entry of `qubits_applied`; for example, `"X,Z"` on `[1, 3]` gives ⟨X₁Z₃⟩.
The result is real-valued.

# Examples
```julia
state = apply(zero_state(2), [Op("H", 1), Op("CX", 1, 2)])
rho = to_rho(state)
@assert correlation(state, "Z,Z", [1, 2]) ≈ 1
@assert correlation(rho, "X,X", [1, 2]) ≈ 1
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
    mag_moments(rho::AbstractMatrixS, op_str::String, moment_order::Int)

Return `tr(rho * M^moment_order)`, where `M` is the sum of the single-qubit
operator `op_str` over all qubits. This is the total magnetization, without
division by the number of qubits. The moment order must be nonnegative;
order zero returns `tr(rho)`.
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
    mag_moments(m::Measurement, moment_order::Int)
    mag_moments(m::Measurement; max_moment=12)

Compute moments of the total measured magnetization `N - 2*number_of_ones`
using the outcome weights in `m`. The keyword form returns orders 1 through
`max_moment`, in the basis recorded by the measurement.
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
    entanglement_entropy(psi::AbstractVectorS; spectrum_bool=nothing)

Compute the entropy across the cut separating the last `floor(N/2)` qubits
from the rest of a normalized pure state, using natural logarithms.

Return `(entropy, spectrum)` by default. Set `spectrum_bool=false` to return
only the entropy, or `true` for the tuple. The spectrum contains `-log(p)`
for the nonzero Schmidt probabilities `p`.
"""
function entanglement_entropy(psi::AbstractVectorS; spectrum_bool::Union{Nothing,Bool}=nothing)
    N = get_N(psi)
    partA = N ÷ 2
    amplitudes = reshape(psi, (2^partA, 2^(N-partA)))
    probabilities = la.svdvals(Matrix(amplitudes)).^2
    probabilities = probabilities[probabilities .> 0]
    entropy = -sum(probabilities .* log.(probabilities))
    return spectrum_bool === false ? entropy : (entropy, -log.(probabilities))
end

"""
    entanglement_entropy(rho::AbstractMatrix; spectrum_bool=nothing)

Return `(entropy, spectrum)` for the reduced density matrix of the last
`floor(N/2)` qubits, using natural logarithms. Set `spectrum_bool=false` for
just the entropy, or `true` for the tuple. For a mixed global state, reduced entropy includes
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
