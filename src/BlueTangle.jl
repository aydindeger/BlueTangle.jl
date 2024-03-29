module BlueTangle

include("all.jl")

export get_N, ro3, fields, sample_outcomes, get_probabilities_from_sample,  expand_multi_op, string_to_matrix, get_corr_from_measurement, get_expect_from_measurement, get_expect, get_corr, apply_op, apply_op_rho, classical_shadow, compile, quantum_circuit, sample, circuit_to_state, circuit_to_rho, show_basis
export QuantumOps,Op,ifOp,Measurement,QuantumChannel,Circuit,Options
export hilbert, hilbert_op, state_vector_create, init_state_create, state_to_measurement
export gate,gates1,gates2,random_ops,random_clifford,Noise1,Noise2,apply_noise,U1,U2,U3,is_kraus_valid,apply_twirl,custom_noise,cnot_amplifier!,op_amplifier!,linear_fit,quadratic_fit,error_mitigate_data
export plot_measurement, plot_circuit, savefigure
export entanglement_entropy,clasical_shadow,mag_moments_from_rho
export trotter_ising

"""
Alias:
```
get_N(state::SparseVector) -> Int
get_N(rho::SparseMatrixCSC) -> Int
```

Calculate the number of qubits (`N`) from a quantum state vector or a density matrix.

- `state`: A sparse quantum state vector.
- `rho`: A sparse density matrix.

Returns `N`, the number of qubits, calculated as the logarithm base 2 of the length of the state vector or the size of the density matrix.
"""
get_N(state::SparseVector)=Int(log(2,length(state)))
get_N(rho::SparseMatrixCSC)=Int(log(2,size(rho,1)))

##
##========== sample outcomes ==========

"""
`sample_outcomes(state::SparseVector, shots::Int) -> Vector`

Sample outcomes from a quantum state vector based on the probability distribution.

- `state`: A sparse quantum state vector.
- `shots`: Number of samples to be taken.

Returns a vector of sampled outcomes.
"""
function sample_outcomes(state::SparseVector, shots)

    N=Int(log(2,length(state)))

    # Convert state to probability distribution
    probs = abs2.(state)
    
    bit_basis=0:2^N-1
    # fock_basis=int2bit.(bit_basis,N)
    # mag_basis=replace.(basis,1=>-1,0=>1)
    # str_basis = join.(map.(string, fock_basis))
    
    sampled_outcomes = StatsBase.sample(bit_basis, StatsBase.Weights(probs), shots)
    
    # Return the sampled outcomes
    return sampled_outcomes
end

"""
`get_probabilities_from_sample(sample_outcomes::Vector, N::Int) -> (Vector, Vector)`

Convert a sample of outcomes into probabilities.

- `sample_outcomes`: A vector of sampled outcomes.
- `N`: Number of qubits.

Returns a tuple of vectors: the first vector contains outcomes, and the second vector contains corresponding probabilities.
"""
function get_probabilities_from_sample(sample_outcomes::Vector, N::Int)
    # Preallocate the frequency array
    freq = zeros(Int, 2^N)

    # Increment frequencies
    for outcome in sample_outcomes
        freq[outcome + 1] += 1  # +1 because Julia arrays are 1-indexed
    end

    total_samples = length(sample_outcomes)

    # Convert frequencies to probabilities
    probabilities = freq ./ total_samples

    # Filter out non-zero probabilities and their corresponding outcomes
    nonzero_indices = findall(x -> x > 0, probabilities)
    outcomes = nonzero_indices .- 1  # Adjusting back to 0-indexed outcomes
    probs = probabilities[nonzero_indices]

    return (outcomes, probs)
end

sample_state(state::SparseVector, shots::Int)=get_probabilities_from_sample(sample_outcomes(state, shots),get_N(state))
sample_state_exact(state::SparseVector)=(0:2^get_N(state)-1,abs2.(state))
sample_rho_exact(rho::SparseMatrixCSC)=(collect(0:2^get_N(rho)-1),real(diag(rho)),get_N(rho))

##========== sample outcomes ==========

##
##========== functions ==========

"""
`expand_multi_op(list_of_operators::String, qubits_applied::Vector, N::Int) -> Matrix`

Expand multiple quantum operators over a specified set of qubits.

- `list_of_operators`: A string representing a list of operators.
- `qubits_applied`: A vector of qubits on which operators are applied.
- `N`: Total number of qubits.

Returns a matrix representing the expanded operators.
"""
function expand_multi_op(list_of_operators::String,qubits_applied::Vector,N::Int)

    ops_str =String.(split(list_of_operators,","))
    result = ["I" for _ in 1:N]

    if length(qubits_applied) != length(ops_str)
        throw("qubit number does not match with operators")
    end
    
    for (op, qubit) in zip(ops_str, qubits_applied)
        result[qubit] = op
    end
    
    str_ops=join(result, ",")

    return string_to_matrix(str_ops)
end

"""
`string_to_matrix(list_of_operators::String) -> Matrix`

Convert a string of comma-separated operators into a matrix representation.

- `list_of_operators`: A string representing a list of operators, e.g.: "Z,Z,I,P(.2)"

Returns a matrix representation of the operators.
"""
function string_to_matrix(list_of_operators::String)
    ops_str=String.(split(list_of_operators,","));
    return foldl(kron,sparse.(gates1.(ops_str)))
end


##========== functions ==========

##
##========== expectation and correlations from measurement ==========

"""
`get_corr_from_measurement(m::Measurement, qubits::Vector) -> Any`

Calculate correlation from a measurement on specified qubits.

- `m`: A `Measurement` object.
- `qubits`: Vector of qubits for measuring correlations, e.g. eg: qubits=[1,3] measure <Z1Z3>

Returns the calculated correlation.
"""
get_corr_from_measurement(m::Measurement,qubits::Vector)=_sample_to_expectation((m.fock_basis,m.sample),qubits)

"""
Alias:
```
get_expect_from_measurement(m::Measurement, qubit::Int) -> Float64
get_expect_from_measurement(m::Measurement) -> Vector
```

Calculate the expectation value from a measurement.

- `m`: A `Measurement` object.
- `qubit`: (Optional) Specific qubit for which to calculate the expectation value.

Returns the expectation value(s).
"""
get_expect_from_measurement(m::Measurement,qubit::Int)=_sample_to_expectation(m::Measurement,[qubit])

get_expect_from_measurement(m::Measurement)=[_sample_to_expectation(m,[i]) for i=1:m.number_of_qubits]

##========== expectation and correlations from measurement ==========

##========== expectation from state & rho and operator ==========

"""
Alias for density matrix:
```
get_expect(state::SparseVector, op::QuantumOps) -> Float64
get_expect(state::SparseVector, op_str::String, qubit::Int) -> Float64
get_expect(state::SparseVector, matrix::SparseMatrixCSC) -> Float64
get_expect(state::SparseVector, op_str::String) -> Vector{Float64}
```

Alias for state vector:
```
get_expect(rho::SparseMatrixCSC, op::QuantumOps) -> Float64
get_expect(rho::SparseMatrixCSC, op_str::String, qubit::Int) -> Float64
get_expect(rho::SparseMatrixCSC, matrix::SparseMatrixCSC) -> Float64
get_expect(rho::SparseMatrixCSC, op_str::String) -> Vector{Float64}
```

Calculate the expectation value for quantum states or density matrices given an operator. This function has several forms depending on the input parameters:

- `get_expect(state::SparseVector, op::QuantumOps)`: Computes the expectation value for a quantum state vector with a given operator.

- `get_expect(rho::SparseMatrixCSC, op::QuantumOps)`: Computes the expectation value for a density matrix with a given operator.

- `get_expect(state::SparseVector, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a quantum state vector with an operator specified as a string.

- `get_expect(rho::SparseMatrixCSC, op_str::String, qubit::Int)`: Computes the expectation value for a specific qubit in a density matrix with an operator specified as a string.

- `get_expect(state::SparseVector, op_str::String)`: Computes the expectation values for all qubits in a quantum state vector given an operator as a string.

- `get_expect(rho::SparseMatrixCSC, op_str::String)`: Computes the expectation values for all qubits in a density matrix given an operator as a string.

- `get_expect(state::SparseVector, matrix::SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a state vector.

- `get_expect(rho::SparseMatrixCSC, matrix::SparseMatrixCSC)`: Computes the expectation value using a sparse matrix representation of an operator for a density matrix.

# Arguments
- `state::SparseVector`: The quantum state vector.
- `rho::SparseMatrixCSC`: The density matrix.
- `op::QuantumOps`: The quantum operator.
- `op_str::String`: The string representation of the operator.
- `qubit::Int`: The specific qubit index.
- `matrix::SparseMatrixCSC`: The sparse matrix representation of the operator.

# Returns
- The expectation value as a `Float64` or a vector of `Float64` for multiple qubits.
"""
get_expect(state::SparseVector,op::QuantumOps)=real(state'*hilbert_op(state,op,get_N(state)))
get_expect(rho::SparseMatrixCSC,op::QuantumOps)=real(tr(rho*hilbert_op(op,get_N(rho))))

get_expect(state::SparseVector,op_str::String,qubit::Int)=real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state)
get_expect(rho::SparseMatrixCSC,op_str::String,qubit::Int)=real(tr(rho*expand_multi_op(op_str,[qubit],get_N(rho))))

get_expect(state::SparseVector,op_str::String)=[real(state'*expand_multi_op(op_str,[qubit],get_N(state))*state) for qubit=1:get_N(state)]
get_expect(rho::SparseMatrixCSC,op_str::String)=[real(tr(rho*expand_multi_op(op_str,[qubit],get_N(rho)))) for qubit=1:get_N(rho)]

get_expect(state::SparseVector,matrix::SparseMatrixCSC)=real(state'*matrix*state)
get_expect(rho::SparseMatrixCSC,matrix::SparseMatrixCSC)=real(tr(rho*matrix))

"""
Alias:
```
get_corr(state::SparseVector, list_of_operators::String, qubits_applied::Vector) -> Float64
get_corr(rho::SparseMatrixCSC, list_of_operators::String, qubits_applied::Vector) -> Float64
```

Calculate the correlation for a given set of operators applied to specific qubits in either a quantum state vector or a density matrix. This function has two primary forms:

- `get_corr(state::SparseVector, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a quantum state vector (`state`) with a specified list of operators and qubits.

- `get_corr(rho::SparseMatrixCSC, list_of_operators::String, qubits_applied::Vector)`: Computes the correlation for a density matrix (`rho`) with a specified list of operators and qubits.

The `corr_from_rho` function is an alias to `get_corr` for density matrices.

# Arguments
- `state::SparseVector`: The quantum state vector.
- `rho::SparseMatrixCSC`: The density matrix.
- `list_of_operators::String`: A string representing a list of operators, e.g., "Z,Z".
- `qubits_applied::Vector`: A vector of qubit indices on which the operators are applied.

# Returns
- `Float64`: The computed correlation value.

# Examples
```julia
# For a state vector
state = SparseVector([...]) # define your state vector
correlation = get_corr(state, "Z,Z", [2, 4])

# For a density matrix
rho = SparseMatrixCSC([...]) # define your density matrix
correlation = get_corr(rho, "Z,Z", [2, 4])
```

"""
function get_corr(state::SparseVector,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(state))
    return real(state'*matrix*state)
end

function get_corr(rho::SparseMatrixCSC,list_of_operators::String,qubits_applied::Vector)
    matrix=expand_multi_op(list_of_operators,qubits_applied,get_N(rho))
    return real(tr(rho*matrix))
end

##========== expectation from state & rho and operator ==========


##***## sample outcomes

##***## apply gate


# function apply_op!(state::SparseVector,op::QuantumOps)
    
#     if op.q!=1 && abs(op.qubit-op.target_qubit)>1
#         throw("non-local gate $(op.name) is not allowed!")
#     end

#     # random measurement mid-circuit
#     if uppercase(op.name)=="MR" || uppercase(op.name)=="M(R)"
#         rOp=rand(1:3)
#         println("random mid-measurement applied in $(["X","Y","Z"][rOp]) basis")
        
#         if typeof(op)==ifOp
#             op=ifOp(["MX","MY","MZ"][rOp],op.qubit,op.if01)
#         else
#             op=Op(["MX","MY","MZ"][rOp],op.qubit)
#         end
#     end

#     # state[:]=_extend_op(op,get_N(state))*state

#     N=get_N(state)
#     res=hilbert_op(state,op,N)
#     state[:]=res

#     if op.noise!=false #if measurement then it applies measurement `channel`
#         _apply_kraus!(state,op)
#     end

# end

"""
`apply_op(state::SparseVector, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `state`: A sparse quantum state vector to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
"""
function apply_op(state::SparseVector,op::QuantumOps)
    
    if op.q!=1 && abs(op.qubit-op.target_qubit)>1
        throw("non-local gate $(op.name) is not allowed!")
    end

    # random measurement mid-circuit
    if uppercase(op.name)=="MR" || uppercase(op.name)=="M(R)"
        rOp=rand(1:3)
        println("random mid-measurement applied in $(["X","Y","Z"][rOp]) basis")
        
        if typeof(op)==ifOp
            op=ifOp(["MX","MY","MZ"][rOp],op.qubit,op.if01)
        else
            op=Op(["MX","MY","MZ"][rOp],op.qubit)
        end
    end

    state=hilbert_op(state,op,get_N(state))

    if op.noise!=false #if measurement then it applies measurement `channel`
        state=_apply_kraus(state,op)
    end

    return state

end


"""
`apply_op_rho(rho::SparseMatrixCSC, op::QuantumOps)`

Apply a quantum gate operation to a state vector in place.

- `rho`: A sparse quantum density matrix to be modified.
- `op`: A `QuantumOps` object representing the gate operation.

Modifies the state vector directly.
"""
function apply_op_rho(rho::SparseMatrixCSC,op::QuantumOps)

    if op.q!=1 && abs(op.qubit-op.target_qubit)>1
        throw("non-local gate $(op.name) is not allowed!")
    end

    # random measurement mid-circuit
    if uppercase(op.name)=="MR" || uppercase(op.name)=="M(R)"
        rOp=rand(1:3)
        println("random mid-measurement applied in $(["X","Y","Z"][rOp]) basis")
        
        if typeof(op)==ifOp
            op=ifOp(["MX","MY","MZ"][rOp],op.qubit,op.if01)
        else
            op=Op(["MX","MY","MZ"][rOp],op.qubit)
        end

    end

    N=get_N(rho)
    e_op=hilbert_op(op,N)

    rho=e_op*rho*e_op'

    if op.noise!=false  #if measurement then it applies measurement `channel`
        rho=_apply_kraus_rho(rho,op)
    end

    return rho
    
end

##***## apply gate

"""
`classical_shadow(circuit::Circuit, number_of_experiment::Int) -> SparseMatrixCSC`

Construct a density matrix (rho) from classical shadow representation.

- `circuit`: A `Circuit` object.
- `number_of_experiment`: Number of experiments to run.

Returns a sparse density matrix representing the classical shadow.
"""
function classical_shadow(circuit::Circuit,number_of_experiment::Int)

    N=circuit.stats.number_of_qubits
    rho_construct=spzeros(ComplexF64,2^N,2^N)
    
    for _ =1:number_of_experiment
        state=circuit_to_state(circuit)

        ## final random measurement
        measurement_basis=Vector{Int}(undef,N)
        for qubit=1:N
            rOp=rand(1:3)
            state=apply_op(state,Op(["MX","MY","MZ"][rOp],qubit))
            measurement_basis[qubit]=rOp
        end

        classical_state=sample_outcomes(state,1)#shot=1 #this should stay as 1 for noisy experiments

        m_list=sparse.([gate.H,gate.HSp,gate.I])
        classical_bit=int2bit(classical_state[1],N)
        rho_single_list=[]
        
        for i in 1:N
            bv = _b2v(classical_bit[i])
            m = measurement_basis[i]
            Um = m_list[m]
            
            # Perform operations more efficiently (check if they can be optimized further)
            push!(rho_single_list, 3 * (Um' * bv * bv' * Um) - sparse(gate.I))
        end

        rho_construct += foldl(kron,rho_single_list)/number_of_experiment

    end

    return rho_construct
end

##***## compile

"""
`quantum_circuit(ops::Vector{T}, options::Options=Options()) -> Circuit where T <: QuantumOps`

Alias for [`compile`](@ref)
"""
quantum_circuit(ops::Vector{T}, options::Options=Options()) where T <: QuantumOps=compile(ops, options)

"""
Alias:
```
quantum_circuit(ops::Vector{T}, options::Options=Options()) -> Circuit where T <: QuantumOps
compile(ops::Vector{T}, options::Options=Options()) -> Circuit where T <: QuantumOps
get_circuit_from_ops(ops::Vector{T}, options::Options=Options()) -> Circuit where T <: QuantumOps
```

Compile a set of quantum operations into a circuit.

- `ops`: A vector of `QuantumOps` objects.
- `options`: Optional compilation options.

Returns a compiled `Circuit` object.
"""
function compile(ops::Vector{T}, options::Options=Options()) where T <: QuantumOps

    local_ops=Vector{QuantumOps}()
    
    two_qubit_count=0
    cx_count=0
    mid_measurement_count=0

    for op in ops

        ## stats
        if op.name=="CX" || op.name=="CNOT"
            cx_count+=1
        end
        
        if op.q==2
            two_qubit_count+=1
        end

        if _is_it_measurement(op.name)
            mid_measurement_count+=1
        end
        ##

        #overall circuit noise
        if (op.q==1 && options.noise1==false) || (op.q==2 && options.noise2==false)

            op_n=op #new_operator

        else

            if op.noise!=false #individual operator has noise so do not override with global options.
                op_n=op
            else
                if op.q==1
                    op_n=Op(op.name,op.mat,op.qubit,options.noise1)
                elseif op.q==2
                    op_n=Op(op.name,op.mat,op.qubit,op.target_qubit,options.noise2)
                else
                    throw("only works for 1 or 2 qubit operators")
                end
            end

        end

        #two qubit - check locality add swaps
        if op_n.q==2 
            distance=abs(op_n.qubit - op_n.target_qubit)
            if distance > 1 #not local
                println("Nonlocal operation warning! Swap will be inserted.")
                append!(local_ops,_non_local_gates(op_n;swap_error=options.swap_error))
            elseif distance==0
                throw("control and target qubits cannot be same!")
            else
                push!(local_ops,op_n) #already local
            end
        elseif op_n.q==1 #single qubit
            push!(local_ops,op_n) #already local
        end

    end

    number_of_qubits=maximum([max(op.qubit, typeof(op)==Op ? op.target_qubit : op.qubit) for op in local_ops])

    len=length(local_ops)
    depth=_calculate_circuit_depth(local_ops)

    stats=(number_of_qubits=number_of_qubits,gate_count=len,depth=depth,cx_count=cx_count,two_qubit_count=two_qubit_count,mid_measurement_count=mid_measurement_count)
    return Circuit(stats,options,local_ops)
end

"""
`sample(circuit::Circuit, number_of_experiment::Int) -> Measurement`

Measure a quantum circuit multiple times.

- `circuit`: A `Circuit` object.
- `number_of_experiment`: Number of times to execute the circuit.

Returns a `Measurement` object with the results.
"""
function sample(circuit::Circuit,number_of_experiment::Int)

    if circuit.options.zne==true
        measurements=Vector{Measurement}()
        for id=0:3 #id=0 means no ZNE
            # println("*** ZNE measurements - extra $(id) CNOT pair applied ***")
            push!(measurements,sample(circuit::Circuit,number_of_experiment,id))
        end
        return measurements
    else
        return sample(circuit::Circuit,number_of_experiment::Int,0)
    end

end

function _measurement_mitigate_inv_matrix(N::Int,final_measurement_error::QuantumChannel)
    mat=[]
    for a=0:2^N-1

        state=spzeros(2^N)
        state[a+1]=1

        rho0=state*state'
        for i=1:N
            op=Op("MZ",gate.I,i,final_measurement_error)
            rho0=apply_op_rho(rho0,op)
        end
        push!(mat,Vector(diag(rho0)))
    end
    return inv(hcat(mat...))
end


"""
`state_to_measurement(state::SparseVector,number_of_experiment::Int)`

this creates a measurement object from state vector.
"""
function state_to_measurement(state::SparseVector,number_of_experiment::Int)

    all_sample=[]
    N=get_N(state)
    rho_construct=spzeros(ComplexF64,2^N,2^N)

    int_basis_unsorted,avg_prob_unsorted=sample_state(state,number_of_experiment)
    sorted_pos=sortperm(int_basis_unsorted) #sort integer basis
    int_basis=int_basis_unsorted[sorted_pos]
    avg_prob=avg_prob_unsorted[sorted_pos]

    fock=int2bit.(int_basis,N)
    expect=[_sample_to_expectation((fock,avg_prob),[i]) for i=1:N]
    mag_moments=[mag_moments_from_measurement(N,int_basis,avg_prob,moment_order) for moment_order=1:12]

    return Measurement(int_basis,fock,avg_prob,expect,mag_moments,"0",number_of_experiment,"state to measurement",N,rho_construct)
    
end



function sample(circuit::Circuit,number_of_experiment::Int,id::Int)

    all_sample=[]
    N=circuit.stats.number_of_qubits

    # measurement_basis=zeros(Int,number_of_experiment)

    rho_construct=spzeros(ComplexF64,2^N,2^N)
    
    for _ =1:number_of_experiment
        state=circuit_to_state(circuit;id=id)
        _final_measurement!(state,circuit.options)#notice how state is changed

        classical_state=sample_outcomes(state,1)#shot=1 #this should stay as 1 for noisy experiments
        append!(all_sample,classical_state)

        if circuit.options.density_matrix==true
            rho_construct += (state*state')/number_of_experiment
        end

    end

    int_basis_unsorted,avg_prob_unsorted=get_probabilities_from_sample(all_sample,N)
    sorted_pos=sortperm(int_basis_unsorted) #sort integer basis
    int_basis=int_basis_unsorted[sorted_pos]
    avg_prob=avg_prob_unsorted[sorted_pos]

    if circuit.options.measurement_mitigate==true && circuit.options.final_measurement_error != false
        println("measurement error mitigation applied")
        avg_prob=_measurement_mitigate_inv_matrix(N,circuit.options.final_measurement_error)*avg_prob
    end

    fock=int2bit.(int_basis,N)
    expect=[_sample_to_expectation((fock,avg_prob),[i]) for i=1:N]
    mag_moments=[mag_moments_from_measurement(N,int_basis,avg_prob,moment_order) for moment_order=1:12]

    return Measurement(int_basis,fock,avg_prob,expect,mag_moments,circuit.options.measurement_basis,number_of_experiment,circuit.options.circuit_name,N,rho_construct)
    
end


"""
`cnot_amplifier!(ops::Vector{T}, CNOT_pair=0) where T <: QuantumOps`

Amplify the presence of CNOT (or CX) operations in a vector of quantum operations.

This function adds extra pair of CNOT operation in the `ops` vector a specific number of times in place. This is useful for amplifying the effect of noise via CNOT operations in a sequence of quantum operations.

# Arguments
- `ops::Vector{T}`: A vector of quantum operations, where `T` is a subtype of `QuantumOps`.
- `CNOT_pair::Int` (optional): The number of additional pairs of CNOT operations to insert for each original CNOT operation in `ops`. The default value is 0, which means no additional operations are inserted.

# Examples
```julia
ops = [Op("H",1), Op("CNOT",1,2), Op("X",2)]
cnot_amplifier!(ops, 1)
# `ops` will be modified to: [Op("H",1), Op("CNOT",1,2), Op("CNOT",1,2), Op("CNOT",1,2), Op("X",2)]
```
"""
function cnot_amplifier!(ops::Vector{T},CNOT_pair=0) where T <: QuantumOps

    i = length(ops)
    while i > 0
        if ops[i].name=="CNOT" || ops[i].name=="CX"
            for _ = 1:2CNOT_pair
                insert!(ops, i, ops[i])
            end
        end
        i -= 1
    end

end

"""
`op_amplifier!(ops::Vector{T},op_pair=0) where T <: QuantumOps`

same as `cnot_amplifier!` but for all operations
"""
function op_amplifier!(ops::Vector{T},op_pair=0) where T <: QuantumOps

    i = length(ops)
    while i > 0
        for _ = 1:2op_pair
            insert!(ops, i, ops[i])
        end
        i -= 1
    end

end



"""
`circuit_to_state(circuit::Circuit; init_state::SparseVector=sparse([])) -> SparseVector`

Convert a quantum circuit to a state vector.

- `circuit`: A `Circuit` object.
- `init_state`: (Optional) Initial state vector.

Returns a state vector representing the circuit.
"""
function circuit_to_state(circuit::Circuit;init_state::SparseVector=sparse([]),id::Int=0)

    N=circuit.stats.number_of_qubits

    if isempty(init_state)
        state=init_state_create(N)
    else
        if N!=get_N(init_state)
            throw("number of qubits error")
        end
        state=copy(init_state)
    end

    for op in circuit.ops

            if op.q==2 && circuit.options.zne==true && (op.name=="CNOT" || op.name=="CX") && id>0 #apply ZNE
            
                for cnot_pair=1:2id+1 #number of id = extra CNOT pair

                    if circuit.options.twirl==true #for each CNOT, twirl applies
                        t_ops = apply_twirl(op)
                        for t_op in t_ops
                            state=apply_op(state,t_op)
                        end
                    else #twirl false but ZNE still applies
                        state=apply_op(state,op)
                    end

                end
    
            elseif op.q==2 && circuit.options.twirl==true #apply twirling/note even ZNE true and id==0
    
                t_ops = apply_twirl(op)
                for t_op in t_ops
                    state=apply_op(state,t_op)
                end

            else #no twirling or ZNE
                state=apply_op(state,op)
            end

    end

    return state

end

"""
`circuit_to_rho(circuit::Circuit) -> SparseMatrixCSC`
`get_rho_from_circuit`
`get_state_from_circuit`

Convert a quantum circuit to a density matrix (rho).

- `circuit`: A `Circuit` object.

Returns a density matrix representing the circuit.
"""
function circuit_to_rho(circuit::Circuit;id::Int=0)

    N=circuit.stats.number_of_qubits
    state=init_state_create(N)
    rho=state*state'

    for op in circuit.ops

        if op.q==2 && circuit.options.zne==true && (op.name=="CNOT" || op.name=="CX") && id>0 #apply ZNE
            
            for cnot_pair=1:2id+1 #number of id = extra CNOT pair

                if circuit.options.twirl==true #for each CNOT, twirl applies
                    t_ops = apply_twirl(op)
                    for t_op in t_ops
                        rho=apply_op_rho(rho,t_op)
                    end
                else #twirl false but ZNE still applies
                    rho=apply_op_rho(rho,op)
                end

            end

        elseif op.q==2 && circuit.options.twirl==true #apply twirling/note even ZNE true and id==0

            t_ops = apply_twirl(op)
            for t_op in t_ops
                rho=apply_op_rho(rho,t_op)
            end

        else #no twirling or ZNE
            rho=apply_op_rho(rho,op)
        end

    end

    return rho
end




end


#todo parse following ops = ["H 1", "CNOT 1 2"]