##
##========== Struct ==========

fields(m)=fieldnames(typeof(m))
attributes(m)=fields(m)

"""
`QuantumChannel(q::Int, model::String, p::Float64)`

Represents a quantum noise channel.

- `q`: Number of qubits affected by the noise channel.
- `model`: The noise model type as a string.
- `p`: Probability parameter for the noise model.
- `kraus`: A vector of matrices representing Kraus operators for the channel.

Constructs a QuantumChannel object for specified qubits, noise model, and probability.
"""
struct QuantumChannel
    q::Int
    model::String
    p::Float64
    kraus::Vector{Matrix}
end


"""
`custom_noise(q::Int, name_of_model::String, kraus::Vector{Matrix}) -> QuantumChannel`

Create a custom noise model represented as a `QuantumChannel` object.

# Arguments
- `q::Int`: The number of qubits affected by the noise model. It should be either 1 for single-qubit noise models or 2 for two-qubit noise models.
- `name_of_model::String`: A name for the custom noise model.
- `kraus::Vector{Matrix}`: A vector of matrices representing the Kraus operators for the noise model. 

# Returns
- `QuantumChannel`: A quantum channel object representing the custom noise model.

# Description
This function constructs a `QuantumChannel` object for a custom noise model based on provided Kraus operators. It validates the Kraus operators to ensure they satisfy the necessary conditions for a quantum channel (e.g., trace preservation). If the Kraus operators are not valid, the function throws an error.

# Example
```julia
kraus_ops = [Matrix([[0.9, 0], [0, 1.0]]), Matrix([[0.0, 0.1], [0.0, 0.0]])]
custom_noise_model = custom_noise(1, "MyNoiseModel", kraus_ops)
```

In this example, `custom_noise_model` is a single-qubit noise model named "MyNoiseModel" with specified Kraus operators.

# Notes
Ensure that the provided Kraus operators comply with the rules for quantum channels. Invalid operators will result in an error being thrown by the function.
"""
custom_noise(q::Int,name_of_model::String,kraus::Vector{Matrix})=QuantumChannel(q,name_of_model,is_kraus_valid(kraus) ? kraus : throw("define valid kraus operators"))

"""
`is_kraus_valid(kraus::Vector{Matrix}) -> Bool`

Determines the validity of a set of Kraus operators.

- `kraus`: A vector of matrices, each representing a Kraus operator.

This function checks if the provided Kraus operators form a valid quantum channel. 
It does so by verifying if the sum of the products of each Kraus operator and its adjoint 
(approximately) equals the identity matrix. Returns `true` if the set is valid, `false` otherwise.
"""
function is_kraus_valid(kraus::Vector{Matrix})
    sumk = sum(k' * k for k in kraus)
    sumk ≈ Matrix(I, size(sumk, 1), size(sumk, 1))
end

"""
`Noise1(model::String, p::Float64)`

Constructors for QuantumChannel objects for 1 qubit.

- `model`: The noise model type as a string.
- `p`: Probability parameter for the noise model.

Returns a QuantumChannel object.
"""
Noise1(model::String,p::Float64)=QuantumChannel(1,model,p)


"""
`Noise2(model::String, p::Float64)`

Constructors for QuantumChannel objects for 2 qubits.

- `model`: The noise model type as a string.
- `p`: Probability parameter for the noise model.

Returns a QuantumChannel object.
"""
Noise2(model::String,p::Float64)=QuantumChannel(2,model,p)

function QuantumChannel(q::Int,model::String,p::Float64)
    model=lowercase(model)
    if q==1
        return QuantumChannel(q,model,p,noise_model1(model,p))
    elseif q==2    
        return QuantumChannel(q,model,p,noise_model2(model,p))
    else
        throw("Noise models are available only for 1 and 2 qubits!")
    end
end

##

"""
`QuantumOps`

Abstract type representing quantum operations.
"""
abstract type QuantumOps end

"""
`Op(q::Int, name::String, gate::Matrix, qubit::Int, target_qubit::Int, noise::QuantumChannel) <: QuantumOps`

Represents a quantum operation.

- `q`: Number of qubits involved in the operation.
- `name`: Name of the operation.
- `gate`: Matrix representation of the quantum operation.
- `qubit`: Index of the target qubit.
- `target_qubit`: Index of the target qubit for two-qubit operations.
- `noise`: Noise model associated with the operation.

Constructs an Op object representing a quantum operation with optional noise.
"""
struct Op <: QuantumOps
    q::Int
    name::String
    gate::Matrix
    qubit::Int
    target_qubit::Int
    noise::Union{QuantumChannel, Bool}
end

"""
`ifOp(q::Int, name::String, gate::Matrix, qubit::Int, if01::Tuple{Matrix,Matrix}, noise::QuantumChannel) <: QuantumOps`

Represents a conditional quantum operation used for mid-circuit born measurements. It is specifically designed for mid-circuit measurements in the X, Y, Z, or a random basis (R). Depending on the measurement outcomes (0 or 1), different gates specified in `if01` can be applied.

# Fields
- `q`: Integer representing the number of qubits involved in the operation.
- `name`: String indicating the name of the operation. Valid names are "MX", "MY", "MZ", or "MR", corresponding to operations in the X, Y, Z basis, or a random basis (R), respectively.
- `gate`: Matrix representing the quantum operation.
- `qubit`: Integer specifying the index of the target qubit.
- `if01`: Tuple of two matrices. The first matrix is applied if the measured state of the qubit is 0, and the second matrix is applied if the measured state is 1.
- `noise`: Represents a 'born measurement quantum channel', indicating the noise model associated with the operation.

# Usage
`ifOp` is constructed to represent quantum operations that are conditional on the measurement outcome of a qubit. This operator is particularly useful in quantum circuits for implementing dynamic responses based on mid-circuit measurement results.

# Examples
```julia
# Example of creating an ifOp for a conditional operation based on the X-basis measurement (MX)
conditional_op = ifOp("MX", qubit, (operation_if_0, operation_if_1))
"""
struct ifOp <: QuantumOps
    q::Int
    name::String
    gate::Matrix
    qubit::Int
    if01::Tuple{Matrix,Matrix}
    noise::QuantumChannel

    ifOp(name::String,qubit::Int,if01::Tuple{Matrix,Matrix}) = new(1, name, gates1(name), qubit, if01, _is_it_measurement(name) ?  Noise1(name,0.0) : throw("select MX or MY or MZ or MR basis."))
    ifOp(name::String,qubit::Int) = ifOp(name,qubit,(gate.I,gate.I))
    ifOp(name::String,qubit::Int,if0::Matrix,if1::Matrix)=ifOp(name,qubit,(if0,if1))

end

function Op(name::String,gate::Matrix,qubit::Int,noise::Union{QuantumChannel, Bool})

    if size(gate,1)!=2
        throw("gate size and selected qubits do not match!")
    end

    if noise==false
        return Op(1,name,gate,qubit,-1,noise)
    elseif noise==true
        throw("provide valid noise model")
    else
        if noise.q==1
            return Op(1,name,gate,qubit,-1,noise)
        else
            throw("gate or noise matrix size error!")
        end
    end

end

function Op(name::String,gate::Matrix,qubit::Int)
    return Op(name,gate,qubit,false)
end

"""
todo define function for namer
"""
function Op(name::String,qubit::Int)
    if _is_it_measurement(name)
        born=Noise1(name,0.0) #note that this is a born measurement `channel`
        return Op(name,gates1(name),qubit,born)
    else
       return Op(name,gates1(name),qubit,false)
    end

end

function Op(name::String,qubit::Int,noise::Union{QuantumChannel,Bool})
    if _is_it_measurement(name)
        println("you can't define error with mid-measurements")
        println("we fixed it no worries!")
        Op(name,qubit)
    else
        Op(name,gates1(name),qubit,noise) #if measurement then ignore error
    end
end

_name_with_phase_bool(name::String)=name ∈ ["P","RX","RY","RZ","CP"]

# Op(namePhase::Vector,qubit::Int)=_name_with_phase_bool(namePhase[1]) ? Op("$(namePhase[1])($(namePhase[2]/pi)π)",gates1(namePhase[1],namePhase[2]),qubit,false) : throw("$(namePhase[1]) does not need a parameter")
# Op(namePhase::Vector,qubit::Int,noise::QuantumChannel)=_name_with_phase_bool(namePhase[1]) ? Op("$(namePhase[1])($(namePhase[2]/pi)π)",gates1(namePhase[1],namePhase[2]),qubit,noise) : throw("$(namePhase[1]) does not need a parameter")


function Op(name::String,gate::Matrix,qubit::Int,target_qubit::Int,noise::Union{QuantumChannel,Bool})

    if size(gate,1)!=4 #two qubit operation
        throw("gate size and selected qubits do not match!")
    elseif qubit==target_qubit
        throw("control and target qubits cannot be same!")
    elseif noise==true
        throw("provide valid noise model")
    else
        return Op(2,name,gate,qubit,target_qubit,noise)
    end

    # # Initialize final_gate with a default value
    # if qubit > target_qubit
    #     final_gate = _swap_control_target(gate)
    # elseif qubit==target_qubit
    #     throw("control and target qubits cannot be same!")
    # else
    #     final_gate = gate
    # end

end

Op(name::String,gate::Matrix,qubit::Int,target_qubit::Int)=Op(name,gate,qubit,target_qubit,false)
Op(name::String,qubit::Int,target_qubit::Int)=Op(name,gates2(name),qubit,target_qubit,false)
Op(name::String,qubit::Int,target_qubit::Int,noise::Union{QuantumChannel,Bool})=Op(name,gates2(name),qubit,target_qubit,noise)

# Op(namePhase::Vector,qubit::Int,target_qubit::Int)=_name_with_phase_bool(namePhase[1]) ? Op("$(namePhase[1])($(namePhase[2]/pi)π)",gates2(namePhase[1],namePhase[2]),qubit,target_qubit,false) : throw("gate does not need a parameter")
# Op(namePhase::Vector,qubit::Int,target_qubit::Int,noise::Union{QuantumChannel,Bool})=_name_with_phase_bool(namePhase[1]) ? Op("$(namePhase[1])($(namePhase[2]/pi)π)",gates2(namePhase[1],namePhase[2]),qubit,target_qubit,noise) : throw("gate does not need a parameter")

##========== Struct ==========

"""
`_swap_control_target(matrix::Matrix) -> Matrix`

Swap the control and target qubits in a two-qubit gate matrix.

- `matrix`: A 4x4 matrix representing a two-qubit operation.

Returns a new matrix with swapped control and target qubits.
"""
function _swap_control_target(matrix::Matrix)
    # Ensure the matrix is 4x4
    if size(matrix) != (4, 4)
        throw("Matrix must be 4x4")
    end

    result = zero(matrix)#zeros(Complex{Float64}, 4, 4)
    
    # Define the permutation that swaps the control and target qubits
    permutation = [1, 3, 2, 4]
    
    for i in 1:4
        for j in 1:4
            result[i, j] = matrix[permutation[i], permutation[j]]
        end
    end
    
    return result
end

"""
`_non_local_gates(op::QuantumOps; swap_error::Bool=false) -> Vector{QuantumOps}`

Generate a sequence of operations to implement a non-local gate.

- `op`: The original non-local quantum operation.
- `swap_error`: Optional parameter to include swap errors.

Returns a vector of QuantumOps to implement the non-local gate.
"""
function _non_local_gates(op::QuantumOps;swap_error=false)

    direction = op.qubit < op.target_qubit ? -1 : 1
    ops_swap=Vector{QuantumOps}()

    if direction==-1
        
        for i=op.target_qubit:direction:op.qubit+2
            swap_error==true ? push!(ops_swap, Op("SWAP", i, i-1, op.noise)) : push!(ops_swap, Op("SWAP", i, i-1))
        end
    
        ops=vcat([ops_swap;Op(op.name,op.qubit,op.qubit+1,op.noise);reverse(ops_swap)])
    
    else
    
        for i=op.target_qubit:direction:op.qubit-2
            swap_error==true ? push!(ops_swap, Op("SWAP", i, i+1, op.noise)) : push!(ops_swap, Op("SWAP", i, i+1)) # [[1 0;0 1]]
        end
    
        ops=vcat([ops_swap;Op(op.name,op.qubit,op.qubit-1,op.noise);reverse(ops_swap)])
    
    end
    
    return ops

end

##
##========== Circuit ==========

"""
`Options(circuit_name::String, measurement_basis::String, measurement_error::QuantumChannel, noise1::QuantumChannel, noise2::QuantumChannel, twirl::Bool, swap_error::Bool, density_matrix::Bool)`

Represents configuration options for a quantum circuit.

- `circuit_name`: Name of the circuit.
- `measurement_basis`: Measurement basis used in the circuit.
- `measurement_error`: Noise model for measurement error.
- `noise1`: Noise model for single-qubit operations.
- `noise2`: Noise model for two-qubit operations.
- `twirl`: Boolean flag for twirling operations.
- `swap_error`: Boolean flag for swap errors.
- `density_matrix`: Boolean flag to indicate if density matrix should be calculated.

Constructs an Options object with specified settings for a quantum circuit.
"""
struct Options
    circuit_name::String
    measurement_basis::String
    measurement_error::Union{QuantumChannel, Bool}
    noise1::Union{QuantumChannel,Bool}
    noise2::Union{QuantumChannel,Bool}
    twirl::Bool
    zne::Bool
    swap_error::Bool
    density_matrix::Bool

    # Constructor with default values
    Options(;
        circuit_name="circuit", 
        measurement_basis="Z",
        measurement_error=false,
        noise1=false, 
        noise2=false,
        twirl=false, 
        zne=false,
        swap_error=false,
        density_matrix=false
    ) = new(circuit_name, uppercase(measurement_basis), measurement_error, noise1, noise2, twirl, zne, swap_error, density_matrix)

    Options()
end

"""
`Circuit`

Represents a quantum circuit.

- `stats`: NamedTuple containing statistics about the circuit.
- `options`: Options object with circuit configurations.
- `ops`: Vector of QuantumOps representing the operations in the circuit.

Constructs a Circuit object representing a quantum circuit.
"""
struct Circuit
    stats::NamedTuple
    options::Options
    ops::Vector{QuantumOps}
end

"""
`Measurement(int_basis::Union{Vector, UnitRange}, fock_basis::Vector{Vector{Int}}, sample::Vector, expect::Vector, mag_moments::Vector, measurement_basis::String, number_of_experiment::Int, circuit_name::String, number_of_qubits::Int, density_matrix::SparseMatrixCSC)`

Represents the result of quantum measurements.

- `int_basis`: Basis of measurement represented as integers or a range.
- `fock_basis`: Fock basis representation of the measurement.
- `sample`: Vector of probabilities.
- `expect`: Expectation values of the measurement.
- `mag_moments`: Magnetic moments obtained from the measurement.
- `measurement_basis`: Measurement basis used.
- `number_of_experiment`: Number of experiments performed.
- `circuit_name`: Name of the circuit used.
- `number_of_qubits`: Number of qubits involved.
- `density_matrix`: Density matrix obtained from the measurement.

Constructs a Measurement object to store the results of quantum measurements.
"""
struct Measurement
    int_basis::Union{Vector, UnitRange}
    fock_basis::Vector{Vector{Int64}}
    sample::Vector #probability
    expect::Vector
    mag_moments::Vector
    measurement_basis::String
    number_of_experiment::Int
    circuit_name::String
    number_of_qubits::Int
    density_matrix::SparseMatrixCSC
end

##========== Circuit ==========

"""
`Data`

Represents a collection of quantum computational data.

- `circuit`: Circuit object representing the quantum circuit.
- `measurement`: Measurement object containing results of the quantum measurements.
- `other`: Vector for storing additional data.

Constructs a Data object for storing and managing quantum computational data.
"""
struct Data
    circuit::Circuit
    measurement::Measurement
    other::Vector
end

"""
`save_data(circuit::Circuit, measurement::Measurement, other::Vector=[])`
`load_data(name::String)`

Functions to save and load quantum computational data.

- `save_data`: Serializes and saves a Data object.
- `load_data`: Deserializes and loads a Data object from a file.

Use these functions for data persistence in quantum computations.
"""
save_data(circuit::Circuit,measurement::Measurement,other=[])=serialize("$(circuit.name)$(rand(1:1000)).ser",Data(circuit,measurement,other))
load_data(name::String)=deserialize("$(name).ser")