"""
`_extend_op(op::QuantumOps, N::Int) -> Matrix`

Extends a quantum operation to the full Hilbert space.

- `op`: The quantum operation (QuantumOps object) to be extended.
- `N`: Total number of qubits in the system.

Returns a matrix representation of the extended operation, ensuring it acts on the entire Hilbert space of `N` qubits.
"""
function _extend_op(op::QuantumOps,N::Int)
    id = [1 0; 0 1]  # Identity

    if op.qubit>N #|| op.target_qubit>N
        throw("qubit cannot be larger than total system size")
    end

    if op.q==1
        e_ops=[x==op.qubit ? sparse(op.gate) : sparse(id) for x=1:N]
    elseif op.q==2
        final_gate=op.qubit > op.target_qubit ? sparse(_swap_control_target(op.gate)) : op.gate
        e_ops=[x==op.qubit ? final_gate : sparse(id) for x=1:N]
        deleteat!(e_ops,op.target_qubit)#target qubit
    end

    return foldl(kron,e_ops)
end

## todo make this inplace function
"""
`_extend_apply(state::SparseVector, op::QuantumOps) -> SparseVector`

Applies an extended quantum operation to a given quantum state.

- `state`: The quantum state vector (SparseVector) to which the operation is applied.
- `op`: The quantum operation (QuantumOps object) to be applied.

Returns the new state vector after applying the extended operation.
"""
_extend_apply(state::SparseVector,op::QuantumOps)=_extend_op(op,get_N(state))*state

"""
`_extend_kraus(N::Int, op::QuantumOps) -> Vector`

Generates extended Kraus operators for a given quantum operation.

- `N`: Total number of qubits in the system.
- `op`: The quantum operation (QuantumOps object) for which to generate Kraus operators.

Returns a vector of Kraus operators extended to the full Hilbert space.
"""
function _extend_kraus(N::Int,op::QuantumOps)

    id = [1 0; 0 1]  # Identity
    
    if op.q==1
        ek_ops=[]
        for kraus in op.noise.kraus
            push!(ek_ops,foldl(kron,[x==op.qubit ? sparse(kraus) : sparse(id) for x=1:N]))
        end
    elseif op.q==2
        ek_ops=[]
        for kraus in op.noise.kraus
            ops_two=[x==op.qubit ? sparse(kraus) : sparse(id) for x=1:N]
            deleteat!(ops_two,op.target_qubit)
            push!(ek_ops,foldl(kron,ops_two))
        end
    else
        throw("selected qubit is more than 2")
    end

    return ek_ops
end

"""
`_weighted_sample(ek_ops::Vector{Any}, probs::Vector{Float64}) -> Tuple`

Selects a Kraus operator from a set based on their associated probabilities.

- `ek_ops`: Vector of extended Kraus operators.
- `probs`: Vector of probabilities corresponding to each Kraus operator.

Returns a tuple containing the index and the selected Kraus operator.
"""
function _weighted_sample(ek_ops::Vector{Any},probs::Vector{Float64})

    rval=rand()

    for (i, cum_weight) in enumerate(cumsum(probs))
        if rval <= cum_weight
            return i,ek_ops[i]
        end
    end
end

"""
`_apply_kraus!(state::SparseVector, op::QuantumOps)`

Applies a Kraus operator to a quantum state vector in place.

- `state`: The quantum state vector (SparseVector) to be modified.
- `op`: The quantum operation (QuantumOps object) containing the Kraus operators.

Modifies the state vector by applying a randomly selected Kraus operator based on the operation's noise model.
"""
function _apply_kraus!(state::SparseVector,op::QuantumOps)

    N=get_N(state)
    ek_ops=_extend_kraus(N,op)

    probs=[abs(state'*ek_op'ek_op*state) for ek_op=ek_ops]
    # println("probs=",probs)

    if sum(probs) ≈ 1

        # kraus=sample(ek_ops, Weights(probs))
        bit,kraus=_weighted_sample(ek_ops,probs)
        normalization=sqrt(abs(state'*kraus'kraus*state))
        
        state[:]=kraus*state/normalization

        if typeof(op)==ifOp #follow-up gate applied
            ifgate=op.if01[bit]
            # println("$(op.name) measurement on qubit $(op.qubit) resulted in |$(bit-1)>\n$(ifgate) applied)\n")
            state[:]=foldl(kron,[x==op.qubit ? sparse(ifgate) : sparse(gate.I) for x=1:N])*state ##extend apply
        end

        # if !(state ≈ state_n)
        #     throw("error with kraus normalization")
        # end

    else
        throw("kraus operator probability error")
    end
    
    # return ek_ops
end

"""
`_apply_kraus_rho!(rho::SparseMatrixCSC, op::QuantumOps)`

Applies a Kraus operator to a density matrix in place.

- `rho`: The density matrix (SparseMatrixCSC) to be modified.
- `op`: The quantum operation (QuantumOps object) containing the Kraus operators.

Modifies the density matrix by applying the Kraus operators, accounting for quantum noise effects.
"""
function _apply_kraus_rho!(rho::SparseMatrixCSC,op::QuantumOps)

    N=get_N(rho)
    ek_ops=_extend_kraus(N,op)

    if typeof(op)==ifOp #follow-up gate applied after mid-measurement
        # println("$(op.name) measurement on qubit $(op.qubit) performed\n")
        after_ops=[foldl(kron,[x==op.qubit ? sparse(op.if01[i]) : sparse(gate.I) for x=1:N]) for i=1:2]
        rho[:]=sum(after_ops[i] * ek_ops[i] * rho * ek_ops[i]' * after_ops[i]' for i=1:2)
    else
        rho[:]=sum(ek * rho * ek' for ek=ek_ops)
    end

end

"""
`_calculate_circuit_depth(ops::Vector{T}) where T <: QuantumOps -> Int`

Calculates the depth of a quantum circuit.

- `ops`: Vector of quantum operations (QuantumOps objects) constituting the circuit.

Returns the depth of the circuit, defined as the maximum number of operations applied to any single qubit.
"""
function _calculate_circuit_depth(ops::Vector{T}) where T <: QuantumOps
    qubit_layers = Dict{Int, Int}()
    circuit_depth = 0

    for op in ops
        max_layer = maximum([get(qubit_layers, qubit, 0) for qubit in op.qubit])
        current_layer = max_layer + 1
        circuit_depth = max(circuit_depth, current_layer)
        for qubit in op.qubit
            qubit_layers[qubit] = current_layer
        end
    end

    return circuit_depth
end

##
##========== expectation from measurement ==========

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

focks=bprob[1]#bit_to(bprob[1],:fock)
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

mag_moments_from_rho(rho::SparseMatrixCSC,op_str::String;max_moment=12)=[mag_moments_from_rho(rho,op_str,ord) for ord=1:max_moment]

"""
`mag_moments_from_rho(rho::SparseMatrixCSC, op_str::String, moment_order::Int) -> Float64`

Calculates the magnetization moments from a given density matrix.

- `rho`: The density matrix (SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `moment_order`: The order of the magnetization moment to compute.

Returns the computed magnetization moment of the specified order.
"""
function mag_moments_from_rho(rho::SparseMatrixCSC,op_str::String,moment_order::Int)
    N=get_N(rho)
    mag_op=sum(_extend_op(Op(op_str,i),N) for i=1:N)
    return real(tr(mag_op^moment_order*rho))
end

"""
`mag_cumulants_from_rho(rho::SparseMatrixCSC, op_str::String, cumulant_order::Int; max_moment::Int=12) -> Float64`

Calculates the magnetization cumulants from a given density matrix.

- `rho`: The density matrix (SparseMatrixCSC) of the quantum system.
- `op_str`: String representation of the operator used for calculating the magnetization.
- `cumulant_order`: The order of the cumulant to compute.
- `max_moment`: (Optional) Maximum order of moments to consider.

Returns the computed magnetization cumulant of the specified order.
"""
function mag_cumulants_from_rho(rho::SparseMatrixCSC,op_str::String,cumulant_order::Int;max_moment=12)
    moments=[mag_moments_from_rho(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments,cumulant_order)
end

function mag_cumulants_from_rho(rho::SparseMatrixCSC,op_str::String;max_moment=12)
    moments=[mag_moments_from_rho(rho,op_str,i) for i=1:max_moment]
    return cumulants_from_moments(moments)
end


"""
calculates average magnetization moments from sample
"""
mag_moments_from_measurement(m::Measurement,moment_order::Int)=mag_moments_from_measurement(m.number_of_qubits,m.int_basis,m.sample,moment_order)
mag_moments_from_measurement(m::Measurement;max_moment=12)=[mag_moments_from_measurement(m::Measurement,ord) for ord=1:max_moment]

function mag_moments_from_measurement(number_of_qubits::Int,int_basis::Union{UnitRange,Vector},sample::Vector,moment_order::Int)
    N=number_of_qubits
    mag=sum.(bit_to(int_basis,N,:mag))
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
    moments=[mag_moments_from_measurement(m,i) for i=1:max_moment]
    return cumulants_from_moments(moments,cumulant_order)
end

"""
`_final_measurement!(state::SparseVector, options::Options)`

Performs the final measurement on a quantum state based on specified options.

- `state`: The quantum state vector (SparseVector) to be measured.
- `options`: Measurement options specifying the basis and error model.

Modifies the state vector to reflect the measurement outcome.
"""
function _final_measurement!(state::SparseVector,options::Options)

    number_of_qubits=get_N(state)

    # measurement_error or random measurement
    if options.measurement_basis=="Z"
        for qubit=1:number_of_qubits
            apply_op!(state,Op("M(Z)",gate.I,qubit,options.measurement_error))
        end
    elseif options.measurement_basis=="X"
        for qubit=1:number_of_qubits
            apply_op!(state,Op("M(X)",gate.H,qubit,options.measurement_error)) #no measurement error
        end
    elseif options.measurement_basis=="Y"
        for qubit=1:number_of_qubits
            apply_op!(state,Op("M(Y)",gate.HSp,qubit,options.measurement_error)) #no measurement error
        end
    elseif options.measurement_basis=="R"
            # random measurement basis
        for qubit=1:number_of_qubits
            rOp=rand(1:3)
            apply_op!(state,Op(["MX","MY","MZ"][rOp],qubit)) #random measurement
        end
    else
        throw("measurement_basis error!")
    end

end


##
##========== Entanglement entropy ==========

"""
`entanglement_entropy(psi::SparseVector) -> Float64`

Calculates the entanglement entropy of a quantum state.

- `psi`: The quantum state vector (SparseVector) for which to calculate the entropy.

Returns the calculated entanglement entropy.
"""
function entanglement_entropy(psi::SparseVector)

    N=get_N(psi)
    
    mat=reshape(psi,(2^(Int(N/2)),2^(Int(N/2))));
    
    spec=svdvals(Matrix(mat)).^2
    
    spec=spec[spec .> 0]
    
    return sum(-spec.*log.(spec)),-log.(spec)
    
end

"""
`partial_trace(rho) -> Matrix`

Computes the partial trace of a density matrix.

- `rho`: The density matrix (SparseMatrixCSC) for which to perform the partial trace.

Returns the resulting matrix after performing the partial trace.

#todo check this! logspec does not seem correct.
"""
function partial_trace(rho)
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
`entanglement_entropy(rho::SparseMatrixCSC) -> Float64`

Calculates the entanglement entropy of a density matrix.

- `rho`: The density matrix (SparseMatrixCSC) for which to calculate the entropy.

Returns the calculated entanglement entropy.
"""
function entanglement_entropy(rho::SparseMatrixCSC)
    mat=partial_trace(rho)#trace half
    spec=svdvals(mat)
    spec=spec[spec .> 0]
    return sum(-spec.*log.(spec)),-log.(spec) 
end

"""
`error_mitigate_data(xdata::Vector, ydata::Vector)`

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
    # xdata=collect(1:2:2length(ydata))
    model(x, p) = p[1] .+ p[2] .* x
    rfit = LsqFit.curve_fit(model, xdata, ydata, [.5, 1.0])
    est=rfit.param[1]
    se=LsqFit.standard_errors(rfit)[1]
    fit_plot=(0:xdata[end],[model(x,rfit.param) for x=0:xdata[end]])
    est,se,fit_plot
end