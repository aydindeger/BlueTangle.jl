"""
    get_N(state::AbstractVectorS) -> Int
    get_N(rho::AbstractMatrixS) -> Int
"""
get_N(state::AbstractVectorS)=Int(log(2,length(state)))
get_N(rho::AbstractMatrixS)=Int(log(2,size(rho,1)))

get_M(state::it.MPS)=it.siteinds(state)
get_N(psi::it.MPS)=length(get_M(psi))

# Scheduling must reserve every qubit a conditional operation can touch.
function _operated_qubits(op::QuantumOps)
    if isa(op, OpQC) && op.q == 3
        return op.qubit:op.qubit+2
    end
    return filter(>(0), (op.qubit, op.target_qubit, op.control))
end

function _operated_qubits(op::OpF)
    if op.data isa Vector{<:QuantumOps}
        return unique([q for nested in op.data for q in _operated_qubits(nested)])
    elseif op.data isa AbstractMatrix
        return 1:get_N(op.data)
    end
    return ()
end

function _operated_qubits(op::ifOp)
    qubits = [op.qubit]
    for branch in op.if01, branch_op in branch
        append!(qubits, _operated_qubits(branch_op))
    end
    return unique!(qubits)
end

"""
qubit number from a set of operations
"""
function get_N(ops::Vector{<:QuantumOps})
    max_index = 0
    for op in ops, qubit in _operated_qubits(op)
        max_index = max(max_index, qubit)
    end
    return max_index
end

"""
    sample(state::AbstractVectorS, shots::Int) -> Vector

Sample outcomes from a quantum state vector based on the probability distribution.

- `state`: A dense or sparse quantum state vector.
- `shots`: Number of samples to be taken.

Returns a vector of sampled outcomes.
"""
function sample(state::AbstractVectorS, shots)

    N=get_N(state)

    # Convert state to probability distribution
    probs = abs2.(state)
    
    bit_basis=0:2^N-1
    
    sampled_outcomes = sb.sample(bit_basis, sb.Weights(probs), shots)
    
    # Return the sampled outcomes
    return sampled_outcomes
end

function sample(state::sa.SparseVector, shots)
    get_N(state) # Check the Hilbert-space dimension, as for dense states.
    indices, amplitudes = sa.findnz(state)
    # Sample only stored amplitudes, rather than scanning all 2^N basis states.
    return sb.sample(indices .- 1, sb.Weights(abs2.(amplitudes)), shots)
end

sample(psi::it.MPS, shots)=bin2int.(sample_bit(psi,shots))

"""
    get_probs_from_sample(sample::Vector, N::Int) -> (Vector, Vector)

Convert a sample of outcomes into probabilities.

- `sample`: A vector of sampled outcomes.
- `N`: Number of qubits.

Returns a tuple of vectors: the first vector contains outcomes, and the second vector contains corresponding probabilities.
"""
function get_probs_from_sample(sample::Vector, N::Int)
    T = N < 8sizeof(Int) ? Int : BigInt
    counts = Dict{T,Int}()
    for outcome in sample
        oi = T(outcome)
        counts[oi] = get(counts, oi, 0) + 1
    end

    total_samples = length(sample)
    outcomes = collect(keys(counts))
    probs = [counts[o] / total_samples for o in outcomes]

    return sort_vectors(outcomes, probs)
end

function sort_vectors(vector1::Vector,vector2::Vector)
    sorted_pos=sortperm(vector1) #sort integer basis
    return (vector1[sorted_pos],vector2[sorted_pos])
end

sample_state(state::Union{AbstractVectorS,AbstractMatrixS}, shots::Int)=get_probs_from_sample(sample(state, shots),get_N(state))
sample_state(psi::it.MPS, shots::Int)=get_probs_from_sample(sample(psi, shots),get_N(psi))

function sample_exact(state::AbstractVectorS)
    probs = abs2.(state)
    indices = findall(!iszero, probs)
    return indices .- 1, probs[indices]
end

function sample_exact(state::sa.SparseVector)
    a,b=sa.findnz(abs2.(state))
    return a .- 1,b
end

"""
    sample_exact(psi::it.MPS)
"""
function sample_exact(psi::it.MPS; max_qubits::Int=16)
    N = get_N(psi)
    N <= max_qubits || throw(ArgumentError("exact MPS probabilities require 2^N outcomes; use measure(psi, shots), or explicitly raise max_qubits"))
    # Contract once rather than repeating a full MPS contraction for every basis state.
    return collect(0:2^N-1), Vector(abs2.(to_state(psi)))
end

function sample_exact(rho::AbstractMatrixS)
    probabilities = real.(la.diag(rho))
    indices = findall(!iszero, probabilities)
    return indices .- 1, probabilities[indices]
end

function sample(rho::AbstractMatrixS, shots)
    outcomes, probabilities = sample_exact(rho)
    return sb.sample(outcomes, sb.Weights(probabilities), shots)
end



"""
`shadow(circuit::Circuit, number_of_experiment::Int) -> sa.SparseMatrixCSC`

Construct a density matrix (rho) from classical shadow representation.

- `circuit`: A `Circuit` object.
- `number_of_experiment`: Number of experiments to run.

Returns a sparse density matrix representing the classical shadow.
"""
function shadow(circuit::Circuit,number_of_experiment::Int)

    N=circuit.stats.N
    rho_construct=sa.spzeros(ComplexF64,2^N,2^N)
    m_list=[gate.H,gate.HSP,gate.I]
    
    for _ =1:number_of_experiment
        state=to_state(circuit)

        ## final random measurement
        measurement_basis=Vector{Int}(undef,N)
        for qubit=1:N
            rOp=rand(1:3)
            state=apply(state,Op(["Xbasis","Ybasis","Zbasis"][rOp],m_list[rOp],qubit))#note that this is not a real measurement just rotation
            measurement_basis[qubit]=rOp
        end

        classical_state=sample(state,1)#shot=1 #this should stay as 1 for noisy experiments

        classical_bit=int2bin(classical_state[1],N)
        rho_single_list=[]
        
        for i in 1:N
            bv = _bin2state(classical_bit[i])
            m = measurement_basis[i]
            Um = sa.sparse(m_list[m])
            
            # Perform operations more efficiently (check if they can be optimized further)
            push!(rho_single_list, 3 * (Um' * bv * bv' * Um) - sa.sparse(gate.I))
        end

        rho_construct += foldl(kron,rho_single_list)/number_of_experiment

    end

    if la.tr(rho_construct)≈1
        return rho_construct
    else
        throw("la.tr(rho)!≈1")
    end
end


"""
get statistics about the operations
"""
function get_stats(ops::Vector{<:QuantumOps})
    max_index = get_N(ops)
    len = 0
    cx_count=0
    swap_count=0
    two_qubit_count=0
    mid_measurement_count=0

    for op in ops
        isa(op, OpF) && continue

        len+=1
        oname=uppercase(op.name)

        if oname=="CX" || oname=="CNOT"
            cx_count+=1
        end

        if oname=="SWAP"
            swap_count+=1
        end

        if op.q==2
            two_qubit_count+=1
        end

        if op.type=="🔬"
            mid_measurement_count+=1
        end

    end

    return (N=max_index,len=len,depth=_calculate_circuit_depth(ops),cx_count=cx_count,swap_count=swap_count,two_qubit_count=two_qubit_count,mid_measurement_count=mid_measurement_count)
end

"""
`_calculate_circuit_depth(ops::Vector{<:QuantumOps}) -> Int`

Calculates the depth of a quantum circuit.

- `ops`: Vector of quantum operations (QuantumOps objects) constituting the circuit.

Returns the number of layers needed while preserving dependencies on all operated
qubits. Conditional branches reserve their combined support, and `OpF` is a barrier.
"""
function _calculate_circuit_depth(ops::Vector{<:QuantumOps})
    qubit_layers = Dict{Int, Int}()
    circuit_depth = 0
    barrier_layer = 0

    for op in ops
        if isa(op, OpF)
            circuit_depth += 1
            barrier_layer = circuit_depth
            continue
        end

        qubits = _operated_qubits(op)
        max_layer = barrier_layer
        for qubit in qubits
            max_layer = max(max_layer, get(qubit_layers, qubit, 0))
        end
        current_layer = max_layer + 1
        circuit_depth = max(circuit_depth, current_layer)
        for qubit in qubits
            qubit_layers[qubit] = current_layer
        end
    end

    return circuit_depth
end

"""
get_layers for a given set of ops
"""
function get_layers(ops::Vector{<:QuantumOps})
    layers = Vector{Vector{QuantumOps}}()
    qubit_layers = Dict{Int, Int}()
    barrier_layer = 0

    for op in ops
        if isa(op, OpF)
            # Arbitrary state functions separate all earlier and later operations.
            push!(layers, [op])
            barrier_layer = length(layers)
            continue
        end

        qubits = _operated_qubits(op)
        op_layer = barrier_layer + 1
        for qubit in qubits
            op_layer = max(op_layer, get(qubit_layers, qubit, 0) + 1)
        end

        # Ensure there is enough space for the new layer
        while length(layers) < op_layer
            push!(layers, Vector{QuantumOps}())
        end
        push!(layers[op_layer], op)

        for qubit in qubits
            qubit_layers[qubit] = op_layer
        end
    end

    return layers
end


"""
optimize_simple(layers::Vector{Vector{QuantumOps}})

Cancel exact identities of matching built-in gates. This optimizes the ideal
circuit; applying a noise model to the shorter circuit naturally changes its
noise exposure. Gates with different `noisy` flags are kept distinct.
"""
function optimize_simple(layers::Vector{Vector{QuantumOps}})
    len = 0 # Initialize len outside the loop
    while len != length(layers) # Continue until no change in number of layers
        len = length(layers)
        for (l, layer) in enumerate(layers)
            i = 1
            while i <= length(layer)
                op = layer[i]
                removed = false

                if op isa Op && op.mat isa AbstractMatrix && op.control == -2 && l < length(layers)
                    if op.q == 1 && op.name in ("X", "Y", "Z", "H") && op.mat == gates(op.name)
                        removed = _remove_duplicate(layers, l, i, 1)
                    elseif op.q == 1 && op.name in ("S", "SX", "XSQRT") && op.mat == gates(op.name) && l + 3 <= length(layers)
                        removed = _remove_duplicate(layers, l, i, 3)
                    elseif op.q == 2 && op.name in ("CX", "CNOT", "CY", "CZ", "ECR", "SWAP", "FSWAP") && op.mat == gates(op.name)
                        removed = _remove_two_qubit_duplicate(layers, l, i)
                    end
                end

                if !removed
                    i += 1
                end
            end
        end

        layers = get_layers(vcat(layers...))
        isempty(layers) && break
    end
    return layers
end

function _same_cancellation_gate(a::Op, b::QuantumOps)
    return b isa Op && b.name == a.name && b.qubit == a.qubit &&
           b.target_qubit == a.target_qubit && b.control == a.control &&
           b.noisy == a.noisy && b.mat isa AbstractMatrix && b.mat == a.mat
end

function _remove_duplicate(layers, l, i, num_layers)
    op = layers[l][i]
    for j in 1:num_layers
        next_layer = get(layers, l+j, nothing)
        if isnothing(next_layer)
            return false
        end
        duplicate_pos = findfirst(x -> _same_cancellation_gate(op, x), next_layer)
        if !isa(duplicate_pos, Number)
            return false
        end
    end

    # If we've found duplicates in all required layers, remove them
    deleteat!(layers[l], i)
    for j in 1:num_layers
        deleteat!(layers[l+j], findfirst(x -> _same_cancellation_gate(op, x), layers[l+j]))
    end
    return true
end

function _remove_two_qubit_duplicate(layers, l, i)
    return _remove_duplicate(layers, l, i, 1)
end

"""
    get_optimized_layers(ops::Vector{<:QuantumOps})
"""
get_optimized_layers(ops::Vector{<:QuantumOps})=optimize_simple(get_layers(ops))

"""
    compile(ops::Vector{<: QuantumOps}, options::Options=Options()) -> Circuit

Compile a set of quantum operations into a circuit.

- `ops`: A vector of `QuantumOps` objects.
- `options`: Optional compilation options.

Returns a compiled `Circuit` object.
"""
function compile(ops::Vector{<:QuantumOps}, options::Options=Options(); layout::Union{Layout,Int}=-1)
#implement layers

    stats_init=get_stats(ops)
    
    N=stats_init.N

    #init layout
    if layout==-1
        layout=Layout(sa.sparse(ones(Int,1, N)))
    elseif N > layout.N #!= N
        throw("The number of physical qubits in the layout does not match the requirements of the operations.")
    end

    #add swaps
    local_ops=Vector{QuantumOps}()
    for op in ops
        append!(local_ops,layout.swap(op))#check locality and add swap if necessary
    end

    stats=get_stats(local_ops)
    layers=get_layers(local_ops)

    return Circuit(stats,stats_init,options,layout,layers)
end

# Independent readout channels have a tensor-product classical response matrix.
function _measurement_mitigate_inv_matrix(N::Int, readout_noise::QuantumChannel)
    readout_noise.q == 1 || throw(ArgumentError("readout noise must be a single-qubit channel"))
    inverse = inv(sum(abs2.(k) for k in readout_noise.kraus))
    return foldl(kron, (inverse for _ in 1:N); init=ones(1, 1))
end

function _mitigate_probabilities(bitstr, probabilities, N, readout_noise)
    readout_noise.q == 1 || throw(ArgumentError("readout noise must be a single-qubit channel"))
    inverse = inv(sum(abs2.(k) for k in readout_noise.kraus))
    corrected = zeros(Float64, 2^N)
    corrected[bitstr .+ 1] = probabilities
    # Apply each 2x2 inverse directly: O(N*2^N) work, no 4^N response matrix.
    for shift in 0:N-1
        stride = 1 << shift
        for block in 0:2stride:length(corrected)-1, offset in 1:stride
            i = block + offset
            j = i + stride
            a, b = corrected[i], corrected[j]
            corrected[i] = inverse[1,1]*a + inverse[1,2]*b
            corrected[j] = inverse[2,1]*a + inverse[2,2]*b
        end
    end
    # Finite-shot inversion gives signed estimates. Taking absolute values
    # biases observables and no longer inverts the classical response.
    corrected ./= sum(corrected)
    return collect(0:length(corrected)-1), corrected
end

function _measurement_rho_placeholder(N::Int)
    return N > 16 ? sa.spzeros(ComplexF64, 0, 0) : sa.spzeros(ComplexF64, 2^N, 2^N)
end

function _measurement_from_probs(bitstr, probabilities, N, shots, basis, label, rho)
    fock = int2bin.(bitstr, N)
    expectations = [_sample_to_expectation((fock, probabilities), [i]) for i in 1:N]
    moments = [mag_moments(N, bitstr, probabilities, order) for order in 1:12]
    return Measurement(bitstr, probabilities, expectations, moments, basis, shots, label, N, rho)
end

"""
    measure(state, shots=-1; label="state to measurement")

Measure a statevector or density matrix. `shots=-1` returns exact probabilities.
For an MPS the default is 1000 samples; request `measure(psi, -1)` explicitly for
exact enumeration (limited to 16 qubits). Large sampled MPS results carry an
empty density-matrix placeholder, never an exponentially large zero matrix.
"""
function measure(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, number_of_experiment::Int=-1; label="state to measurement")
    number_of_experiment == -1 || number_of_experiment > 0 || throw(ArgumentError("shots must be positive, or -1 for exact probabilities"))
    N = get_N(state)
    bitstr, probabilities = number_of_experiment == -1 ? sample_exact(state) : sample_state(state, number_of_experiment)
    rho = state isa AbstractMatrix ? sa.sparse(state) : _measurement_rho_placeholder(N)
    return _measurement_from_probs(bitstr, probabilities, N, number_of_experiment, "0", label, rho)
end

measure(psi::it.MPS; label="state to measurement") = measure(psi, 1000; label=label)

function measure(samples::Vector{<:Integer}, N::Int)
    isempty(samples) && throw(ArgumentError("at least one sample is required"))
    bitstr, probabilities = get_probs_from_sample(samples, N)
    return _measurement_from_probs(bitstr, probabilities, N, length(samples), "0", "sample to measurement", _measurement_rho_placeholder(N))
end

_is_cnot_like(op::QuantumOps) = op isa Op && op.q == 2 && op.control == -2 &&
    uppercase(op.name) in ("CNOT", "CX") && op.mat isa AbstractMatrix && op.mat == gate.CX

function _samples_to_measurement(samples::Vector{<:Integer}, circuit::Circuit; rho=_measurement_rho_placeholder(circuit.stats.N))
    N = circuit.stats.N
    bitstr, probabilities = get_probs_from_sample(samples, N)
    options = circuit.options
    if options.measurement_mitigate && options.readout_noise isa QuantumChannel
        bitstr, probabilities = _mitigate_probabilities(bitstr, probabilities, N, options.readout_noise)
    end
    return _measurement_from_probs(bitstr, probabilities, N, length(samples), options.measurement_basis, options.circuit_name, rho)
end

function _circuit_ops(circuit::Circuit)
    operations = QuantumOps[]
    for layer in circuit.layers, op in layer
        if circuit.options.twirl && _is_cnot_like(op)
            append!(operations, apply_twirl(op))
        else
            push!(operations, op)
        end
    end
    return operations
end

function _final_measurement(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, options::Options; kwargs...)
    basis = options.measurement_basis
    basis in ("X", "Y", "Z", "R") || throw(ArgumentError("measurement_basis must be X, Y, Z, or R"))
    for qubit in 1:get_N(state)
        name = basis == "R" ? rand(("H", "HSP", "I")) : basis == "X" ? "H" : basis == "Y" ? "HSP" : "I"
        name != "I" && (state = apply(state, Op(name, qubit); kwargs...))
        if options.readout_noise isa QuantumChannel
            channel = options.readout_noise
            channel.q == 1 || throw(ArgumentError("readout noise must be a single-qubit channel"))
            state = apply(state, OpQC(channel.name, channel.kraus, qubit); kwargs...)
        end
    end
    return state
end

_has_state_callback(op::QuantumOps) = false
_has_state_callback(op::OpF) = op.data isa Function ||
    (op.data isa Vector{<:QuantumOps} && any(_has_state_callback, op.data))

"""
    measure(circuit, shots; backend=:statevector, kwargs...)

Sample final outcomes with statevector, MPS (`:tensor`), or density-matrix
(`:density_matrix`) evolution. MPS accepts `cutoff` and `maxdim`. `run` instead
returns the mid-circuit measurement records.
"""
function measure(circuit::Circuit, number_of_experiment::Int;
                 backend::Union{Symbol,AbstractString}=:statevector, kwargs...)
    number_of_experiment > 0 || throw(ArgumentError("shots must be positive"))
    N = circuit.stats.N
    backend_name = _resolve_backend(backend)
    _unsupported_keyword_check(kwargs, backend_name == "tensor" ? [:cutoff, :maxdim] : Symbol[])
    T = N < 8sizeof(Int) ? Int : BigInt
    all_sample = T[]
    options = circuit.options
    rho = _measurement_rho_placeholder(N)
    if options.density_matrix
        N <= 16 || throw(ArgumentError("density_matrix=true explicitly materializes 4^N entries; disable it for large MPS circuits"))
    end

    ops = [op for layer in circuit.layers for op in layer]
    stochastic = options.noise != false || options.readout_noise != false ||
                 any(op -> op isa OpQC || op isa ifOp || op isa OpF || op.type == "🔬", ops)
    # Density evolution already sums noise and measurement branches exactly.
    reuse_state = (backend_name == "density_matrix" || !stochastic) &&
                  !any(_has_state_callback, ops) && !options.twirl && options.measurement_basis != "R"
    repetitions = reuse_state ? 1 : number_of_experiment
    for _ in 1:repetitions
        state = if backend_name == "tensor"
            apply(_circuit_ops(circuit), zero_state(N_MPS(N)); noise=options.noise, kwargs...)
        elseif backend_name == "density_matrix"
            to_rho(circuit)
        else
            to_state(circuit)
        end
        state = _final_measurement(state, options; kwargs...)
        append!(all_sample, sample(state, reuse_state ? number_of_experiment : 1))
        if options.density_matrix
            vector = state isa it.MPS ? to_state(state) : state
            rho += (vector isa AbstractMatrix ? vector : vector * vector') / repetitions
        end
    end
    return _samples_to_measurement(all_sample, circuit; rho=sa.sparse(rho))
end

function _resolve_backend(backend::Union{Symbol,AbstractString})
    backend_s = lowercase(String(backend))
    if backend_s in ["statevector","state_vector","state","sv"]
        return "statevector"
    elseif backend_s in ["density_matrix","density","dm","rho"]
        return "density_matrix"
    elseif backend_s in ["tensor","tensor_network","tn","mps"]
        return "tensor"
    else
        throw(ArgumentError("Unsupported backend `$backend`. Use `:statevector`, `:tensor`, or `:density_matrix`."))
    end
end

function run(circuit::Circuit,number_of_experiment::Int=1;
             backend::Union{Symbol,AbstractString}=:statevector,kwargs...)
    number_of_experiment > 0 || throw(ArgumentError("`number_of_experiment` must be positive."))

    N=circuit.stats.N
    backend_name = _resolve_backend(backend)
    measurements = Vector{Vector{Int}}(undef,number_of_experiment)

    if backend_name in ("statevector", "density_matrix")
        _unsupported_keyword_check(kwargs,[])
        nm=circuit.options.noise

        for shot=1:number_of_experiment
            state=zero_state(N)
            backend_name == "density_matrix" && (state = state * state')
            ops_shot=_circuit_ops(circuit)
            _,mid_values = apply(ops_shot,state;noise=nm,track_measurements=true)
            measurements[shot] = mid_values
        end
    else
        _unsupported_keyword_check(kwargs,[:cutoff,:maxdim])

        for shot=1:number_of_experiment
            psi=zero_state(N_MPS(N))
            ops_shot=_circuit_ops(circuit)
            _,mid_values = apply(ops_shot,psi;noise=circuit.options.noise,track_measurements=true,kwargs...)
            measurements[shot] = mid_values
        end
    end

    return measurements
end

function _qasm_to_source(qasm_or_path::AbstractString; is_file::Union{Nothing,Bool}=nothing)
    if is_file===true
        return read(qasm_or_path, String)
    elseif is_file===false
        return String(qasm_or_path)
    end

    if isfile(qasm_or_path)
        return read(qasm_or_path, String)
    end

    qasm_text = String(qasm_or_path)
    if occursin(';', qasm_text) || occursin('\n', qasm_text) || occursin("openqasm", lowercase(qasm_text))
        return qasm_text
    end

    throw(ArgumentError("Input is neither an existing file nor recognizable OpenQASM text: `$qasm_or_path`"))
end

function run(ops::Vector{<:QuantumOps},number_of_experiment::Int=1;
             options::Options=Options(),layout::Union{Layout,Int}=-1,
             backend::Union{Symbol,AbstractString}=:statevector,kwargs...)
    circuit = compile(ops, options; layout=layout)
    return run(circuit,number_of_experiment;backend=backend,kwargs...)
end

function run(qasm_or_path::AbstractString,number_of_experiment::Int=1;
             options::Options=Options(),layout::Union{Layout,Int}=-1,
             backend::Union{Symbol,AbstractString}=:statevector,
             is_file::Union{Nothing,Bool}=nothing,kwargs...)
    qasm = _qasm_to_source(qasm_or_path; is_file=is_file)
    ops = from_qasm(qasm)
    return run(ops,number_of_experiment;options=options,layout=layout,backend=backend,kwargs...)
end


"""
`cnot_amplifier!(ops::Vector{<:QuantumOps}, CNOT_pair=0)`

Amplify the presence of CNOT (or CX) operations in a vector of quantum operations.

This function adds extra pair of CNOT operation in the `ops` vector a specific number of times in place. This is useful for amplifying the effect of noise via CNOT operations in a sequence of quantum operations.

# Arguments
- `ops::Vector{<:QuantumOps}`: A vector of quantum operations, where `T` is a subtype of `QuantumOps`.
- `CNOT_pair::Int` (optional): The number of additional pairs of CNOT operations to insert for each original CNOT operation in `ops`. The default value is 0, which means no additional operations are inserted.

# Examples
```julia
ops = [Op("H",1), Op("CNOT",1,2), Op("X",2)]
cnot_amplifier!(ops, 1)
# `ops` will be modified to: [Op("H",1), Op("CNOT",1,2), Op("CNOT",1,2), Op("CNOT",1,2), Op("X",2)]
```
"""
function cnot_amplifier!(ops::Vector{<:QuantumOps},CNOT_pair=0)

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
`op_amplifier!(ops::Vector{<:QuantumOps},op_pair=0)`

same as `cnot_amplifier!` but for all operations
"""
function op_amplifier!(ops::Vector{<:QuantumOps},op_pair=0)

    i = length(ops)
    while i > 0
        for _ = 1:2op_pair
            insert!(ops, i, ops[i])
        end
        i -= 1
    end

end


"""
`to_state(circuit::Circuit; init_state::AbstractVectorS=sa.sparse([])) -> sa.SparseVector`

Convert a quantum circuit to a state vector.

- `circuit`: A `Circuit` object.
- `init_state`: (Optional) Initial state vector.

Returns a state vector representing the circuit.
"""
function to_state(circuit::Circuit; init_state::AbstractVectorS=zero_state(circuit.stats.N))
    get_N(init_state) == circuit.stats.N || throw(DimensionMismatch("initial state must match the circuit qubit count"))
    state=init_state
    ops_shot=_circuit_ops(circuit)
    return apply(ops_shot,state;noise=circuit.options.noise)

end

"""
    to_rho(circuit::Circuit) -> sa.SparseMatrixCSC

Convert a quantum circuit to a density matrix (rho).

- `circuit`: A `Circuit` object.

Returns a density matrix representing the circuit.
"""
function to_rho(circuit::Circuit)
    state = zero_state(circuit.stats.N)
    return apply(_circuit_ops(circuit), state * state'; noise=circuit.options.noise)
end

"""Expand an operation on an N-qubit Hilbert space (legacy alias)."""
hilbert_op(op::QuantumOps, N::Int) = op.expand(N)

"""Compile operations into a circuit (alias for `compile`)."""
quantum_circuit(args...; kwargs...) = compile(args...; kwargs...)


"""
`expand_multi_op(list_of_operators::String, qubits_applied::Vector, N::Int) -> Matrix`

Expand multiple quantum operators over a specified set of qubits.

- `list_of_operators`: A string representing a list of operators.
- `qubits_applied`: A vector of qubits on which operators are applied.
- `N`: Total number of qubits.

Returns a matrix representing the expanded operators.
"""
function expand_multi_op(list_of_operators::String,qubits_applied::Vector{Int},N::Int)

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
`get_operator(list_of_operators::String, qubits_applied::Vector, N::Int) -> Matrix`

Expand multiple quantum operators over a specified set of qubits.

- `list_of_operators`: A string representing a list of operators.
- `qubits_applied`: A vector of qubits on which operators are applied.
- `N`: Total number of qubits.

Returns a matrix representing the expanded operators.
"""
function get_operator(list_of_operators::Union{String,Vector{String}},qubits_applied::Vector{Int},N::Int)
    if isa(list_of_operators, Vector{String})
        list_of_operators=join(list_of_operators, ",")
    end

    return expand_multi_op(list_of_operators,qubits_applied,N)
end

function get_operator(list_of_operators::Union{String,Vector{String}},qubits_applied::Int,N::Int)
    if isa(list_of_operators, Vector{String})
        list_of_operators=join(list_of_operators, ",")
    end
 
    return expand_multi_op(list_of_operators,[qubits_applied],N)
end

"""
    get_operator(observable::Union{String,Vector{String}}, n::Int) -> AbstractMatrix

Constructs the Kronecker product of the operator `observable` repeated `n` times.

# Arguments
- `observable::String`: The operator to be repeated, represented as a string (e.g., "X").
- `n::Int`: The number of times the operator is repeated.

# Returns
- An `AbstractMatrix` representing the Kronecker product of the operator repeated `n` times.

Example: get_operator("X", 3) -> kron(X, X, X)
"""
function get_operator(observable::Union{String,Vector{String}}, N::Int)

    if isa(observable, Vector{String})
        observable=join(observable, ",")
    end

    ops = join(fill(observable, N), ",")
    return string_to_matrix(ops)
end

function get_operator(observable::Union{String,Vector{String}})
    if isa(observable, Vector{String})
        observable=join(observable, ",")
    end
    return string_to_matrix(observable)
end

expand_multi_op(list_of_operators::String,qubits_applied::Int,N::Int)=expand_multi_op(list_of_operators,[qubits_applied],N)

"""
`string_to_matrix(list_of_operators::String) -> Matrix`

Convert a string of comma-separated operators into a matrix representation.

- `list_of_operators`: A string representing a list of operators, e.g.: "Z,Z,sa.IP(.2)"

Returns a matrix representation of the operators.
"""
function string_to_matrix(list_of_operators::String)
    ops_str=String.(split(list_of_operators,","));
    return foldl(kron,sa.sparse.(gates.(ops_str)))
end

"""
    hamming_distance(v1::Vector{Int}, v2::Vector{Int})
"""
hamming_distance(v1::Vector{Int}, v2::Vector{Int})=sum(v1 .⊻ v2)
