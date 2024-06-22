"""
    get_N(state::AbstractVectorS) -> Int
    get_N(rho::sa.SparseMatrixCSC) -> Int
"""
get_N(state::AbstractVectorS)=Int(log(2,length(state)))
get_N(rho::sa.SparseMatrixCSC)=Int(log(2,size(rho,1)))

get_M(state::it.MPS)=it.siteinds(state)
get_N(psi::it.MPS)=length(get_M(psi))

"""
    sample(state::AbstractVectorS, shots::Int) -> Vector

Sample outcomes from a quantum state vector based on the probability distribution.

- `state`: A sparse quantum state vector.
- `shots`: Number of samples to be taken.

Returns a vector of sampled outcomes.
"""
function sample(state::AbstractVectorS, shots)

    N=get_N(state)

    # Convert state to probability distribution
    probs = abs2.(state)
    
    bit_basis=0:2^N-1
    # fock_basis=int2bin.(bit_basis,N)
    # mag_basis=replace.(basis,1=>-1,0=>1)
    # str_basis = join.(map.(string, fock_basis))
    
    sampled_outcomes = sb.sample(bit_basis, sb.Weights(probs), shots)
    
    # Return the sampled outcomes
    return sampled_outcomes
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
    # Preallocate the frequency array
    freq = zeros(Int, 2^N)

    # Increment frequencies
    for outcome in sample
        freq[outcome + 1] += 1  # +1 because Julia arrays are 1-indexed
    end

    total_samples = length(sample)

    # Convert frequencies to probabilities
    probabilities = freq ./ total_samples

    # Filter out non-zero probabilities and their corresponding outcomes
    nonzero_indices = findall(x -> x > 0, probabilities)
    outcomes = nonzero_indices .- 1  # Adjusting back to 0-indexed outcomes
    probs = probabilities[nonzero_indices]

    return sort_vectors(outcomes, probs)
end

function sort_vectors(vector1::Vector,vector2::Vector)
    sorted_pos=sortperm(vector1) #sort integer basis
    return (vector1[sorted_pos],vector2[sorted_pos])
end

sample_state(state::AbstractVectorS, shots::Int)=get_probs_from_sample(sample(state, shots),get_N(state))
sample_state(psi::it.MPS, shots::Int)=get_probs_from_sample(sample(psi, shots),get_N(psi))

function sample_exact(state::AbstractVectorS)
    a,b=sa.findnz(abs2.(state))
    return a .- 1,b
end

"""
    sample_exact(psi::it.MPS)
"""
function sample_exact(psi::it.MPS)

    M=it.siteinds(psi)
    N=length(M)

    probs=[]
    bitstr=collect(0:2^N-1)
    for a=bitstr

        config=BlueTangle.int2bin(a,N) .+ 1
  
        V = it.ITensor(1.)
        for j=1:N
        V *= (psi[j]*it.state(M[j],config[j]))
        end
  
        push!(probs,abs2(it.scalar(V)))
    end

    return bitstr,probs
end


function sample_exact(rho::sa.SparseMatrixCSC)
    a,b=sa.findnz(real(sa.diag(rho)))
    return a .- 1,b
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
    m_list=[gate.H,gate.HSp,gate.I]
    
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
    max_index = 0
    len = 0
    cx_count=0
    swap_count=0
    two_qubit_count=0
    mid_measurement_count=0

    ops=filter(op -> !isa(op, OpF), ops)

    for op in ops

        len+=1
        oname=uppercase(op.name)

        if isa(op, ifOp)
            current_max = op.qubit
        elseif isa(op, OpQC)
            current_max = max(op.qubit, op.target_qubit)
        else
            current_max = max(op.qubit, op.target_qubit, op.control)
        end

        max_index = max(max_index, current_max)

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

Returns the depth of the circuit, defined as the maximum number of operations applied to any single qubit.
"""
function _calculate_circuit_depth(ops::Vector{<:QuantumOps})
    qubit_layers = Dict{Int, Int}()
    circuit_depth = 0

    for op in ops

        # Initialize max_layer based on the control qubit
        max_layer = get(qubit_layers, op.qubit, 0)

        # If there's no -1, update max_layer based on it
        if op.q==2
            max_layer = max(max_layer, get(qubit_layers, op.target_qubit, 0))
        end

        # The current operation will be placed in the next layer
        current_layer = max_layer + 1
        circuit_depth = max(circuit_depth, current_layer)

        # Update the layers for both qubits involved in the operation
        qubit_layers[op.qubit] = current_layer
        if op.q == 2
            qubit_layers[op.target_qubit] = current_layer
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

    # Determine the maximum qubit index
    max_q = maximum([max(o.qubit, _target_find(o), _control_find(o)) for o in ops if !isa(o, OpF)])

    for op in ops
        if isa(op, OpF)
            # Create a new layer for OpF
            push!(layers, [op])
            
            # Reset all qubit layers to ensure separation
            for q in 1:max_q
                qubit_layers[q] = length(layers)
            end

            continue
        end

        # Determine the first possible layer for the operation
        op_layer = 0
        if op.q == 2
            op_layer = max(get(qubit_layers, op.qubit, 0), get(qubit_layers, op.target_qubit, 0), get(qubit_layers, op.control, 0)) + 1
        elseif op.q == 1
            if isa(op, ifOp)
                op_layer = get(qubit_layers, op.qubit, 0) + 1
            else
                op_layer = max(get(qubit_layers, op.qubit, 0), get(qubit_layers, op.control, 0)) + 1
            end
        end

        # Ensure there is enough space for the new layer
        while length(layers) < op_layer
            push!(layers, Vector{QuantumOps}())
        end
        push!(layers[op_layer], op)

        # Update the engagement layer for the qubits involved
        qubit_layers[op.qubit] = op_layer

        if op.q == 2
            qubit_layers[op.target_qubit] = op_layer
        end

        if isa(op, Op) && op.control != -2
            qubit_layers[op.control] = op_layer
        end
    end

    return layers
end


"""
    delete_duplicates(layers::Vector{Vector{QuantumOps}})
"""
function delete_duplicates(layers::Vector{Vector{QuantumOps}})
    
    ##removes obvious repetitions
    for (l,layer)=enumerate(layers)
        for op=layer
            if (op.name=="I" || op.name=="X" || op.name=="Y" || op.name=="Z") && (l<length(layers))
    
                duplicate_pauli_pos=findfirst((x.name==op.name && x.qubit==op.qubit) for x=layers[l+1])
            
                if isa(duplicate_pauli_pos,Number)
                    popat!(layers[l+1],duplicate_pauli_pos)
                end
    
            end
        end
    end

    return get_layers(vcat(layers...))

end

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

"""
`measure_ZNE(circuit::Circuit, number_of_experiment::Int) -> Measurement`

Measure a quantum circuit multiple times.

- `circuit`: A `Circuit` object.
- `number_of_experiment`: Number of times to execute the circuit.

Returns a `Measurement` object with the results.
"""
function measure_ZNE(circuit::Circuit,number_of_experiment::Int,number_of_ZNE::Int=3)

    if isa(circuit.noise,NoiseModel) #if there is any noise then do this
        measurements=Vector{Measurement}()
        for id=0:number_of_ZNE #id=0 means no ZNE
            # println("*** ZNE measurements - extra $(id) CNOT pair applied ***")
            push!(measurements,measure(circuit::Circuit,number_of_experiment,id))
        end
        return measurements
    else
        throw("There is no noise, why ZNE?")
    end

end

function _measurement_mitigate_inv_matrix(N::Int,readout_noise::QuantumChannel)
    vecs=[]
    for a=0:2^N-1

        state=sa.spzeros(2^N)
        state[a+1]=1
        rho0=state*state'

        for qubit=1:N
            rho0=readout_noise.apply(rho0,qubit)
        end

        push!(vecs,Vector(sa.diag(rho0)))
    end
    
    return inv(hcat(vecs...))
end


"""
`measure(state::AbstractVectorS,number_of_experiment::Int)`

this creates a measurement object from state vector.
"""
function measure(state::Union{AbstractVectorS,it.MPS},number_of_experiment::Int=-1;label="state to measurement")

    N=get_N(state)
    rho_construct=sa.spzeros(ComplexF64,2^N,2^N)

    if number_of_experiment==-1
        bitstr,avg_prob=sample_exact(state)
    else
        bitstr,avg_prob=sample_state(state,number_of_experiment)
    end

    fock=int2bin.(bitstr,N)
    expect=[_sample_to_expectation((fock,avg_prob),[i]) for i=1:N]
    mag_moments_list=[mag_moments(N,bitstr,avg_prob,moment_order) for moment_order=1:12]

    return Measurement(bitstr,avg_prob,expect,mag_moments_list,"0",number_of_experiment,label,N,rho_construct)
    
end


function measure(sample::Vector{Int},N::Int)

    rho_construct=sa.spzeros(ComplexF64,2^N,2^N)

    bitstr,avg_prob=get_probs_from_sample(sample,N)

    fock=int2bin.(bitstr,N)
    expect=[BlueTangle._sample_to_expectation((fock,avg_prob),[i]) for i=1:N]
    mag_moments_list=[mag_moments(N,bitstr,avg_prob,moment_order) for moment_order=1:12]

    return Measurement(bitstr,avg_prob,expect,mag_moments_list,"0",length(sample),"sample to measurement",N,rho_construct)
    
end

function measure(circuit::Circuit,number_of_experiment::Int,id::Int=0)

    all_sample=[]
    N=circuit.stats.N

    rho_construct=sa.spzeros(ComplexF64,2^N,2^N)
    
    if circuit.options.noise==false && circuit.stats.mid_measurement_count==0 && circuit.options.readout_noise==false

        state=to_state(circuit,0)# id=0 no noise no zne
        state=_final_measurement(state,circuit.options)#notice how state is changed
        classical_state=sample(state,number_of_experiment) #shots=number_of_experiment
        append!(all_sample,classical_state)
        if circuit.options.density_matrix==true
            rho_construct = state*state'
        end

    else #noisy circuit

        println("noisy or measurement operations detected!")
        for _ =1:number_of_experiment
            state=to_state(circuit,id)
            state=_final_measurement(state,circuit.options)#notice how state is changed

            classical_state=sample(state,1)#shot=1 #this should stay as 1 for noisy experiments
            append!(all_sample,classical_state)

            if circuit.options.density_matrix==true
                rho_construct += (state*state')/number_of_experiment
            end

        end
    end

    bitstr,avg_prob=get_probs_from_sample(all_sample,N)

    if circuit.options.measurement_mitigate==true && isa(circuit.options.readout_noise,QuantumChannel)
        println("measurement error mitigation applied")
        avg_prob=la.normalize(abs.(_measurement_mitigate_inv_matrix(N,circuit.options.readout_noise)*avg_prob),1) #fix size issue here
    end

    fock=int2bin.(bitstr,N)
    expect=[_sample_to_expectation((fock,avg_prob),[i]) for i=1:N]
    mag_moments_list=[mag_moments(N,bitstr,avg_prob,moment_order) for moment_order=1:12] #fix this

    return Measurement(bitstr,avg_prob,expect,mag_moments_list,circuit.options.measurement_basis,number_of_experiment,circuit.options.circuit_name,N,rho_construct)
    
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
function to_state(circuit::Circuit,id::Int=0)

    N=circuit.stats.N
    state=zero_state(N)
    nm=circuit.options.noise #noisemodel

    # for layer in circuit.layers

        # if circuit.options.zne==false && circuit.options.twirl==false #fix this
        #     state=hilbert_layer(N,layer,state)
        # else_final_measurement

            for op=vcat(circuit.layers...)

                if isa(op,OpF)
                    state=apply(state,op;noise=nm)
                    continue
                end

                if id>0 && op.q==2 && (uppercase(op.name)=="CNOT" || uppercase(op.name)=="CX") #apply ZNE
                
                    for cnot_pair=1:2id+1 #number of id = extra CNOT pair

                        if circuit.options.twirl==true #for each CNOT, twirl applies
                            t_ops = apply_twirl(op)
                            for t_op in t_ops
                                state=apply(state,t_op;noise=nm)
                            end
                        else #twirl false but ZNE still applies
                            state=apply(state,op;noise=nm)
                        end

                    end
        
                elseif op.q==2 && circuit.options.twirl==true #apply twirling when id==0
        
                    t_ops = apply_twirl(op)
                    for t_op in t_ops
                        state=apply(state,t_op;noise=nm)
                    end

                else #no twirling or ZNE
                    state=apply(state,op;noise=nm)
                end
            end

        # end
        
    # end

    return state

end

"""
    to_rho(circuit::Circuit) -> sa.SparseMatrixCSC

Convert a quantum circuit to a density matrix (rho).

- `circuit`: A `Circuit` object.

Returns a density matrix representing the circuit.
"""
function to_rho(circuit::Circuit;id::Int=0)

    N=circuit.stats.N
    state=zero_state(N)
    rho=state*state'
    nm=circuit.options.noise #noisemodel

    for op in vcat(circuit.layers...)
        # for op=layer
            if op.q==2 && (op.name=="CNOT" || op.name=="CX") && id>0 #apply ZNE
                
                for cnot_pair=1:2id+1 #number of id = extra CNOT pair

                    if circuit.options.twirl==true #for each CNOT, twirl applies
                        t_ops = apply_twirl(op)
                        for t_op in t_ops
                            rho=apply(rho,t_op;noise=nm)
                        end
                    else #twirl false but ZNE still applies
                        rho=apply(rho,op;noise=nm)
                    end

                end

            elseif op.q==2 && circuit.options.twirl==true #apply twirling/note even ZNE true and id==0

                t_ops = apply_twirl(op)
                for t_op in t_ops
                    rho=apply(rho,t_op;noise=nm)
                end

            else #no twirling or ZNE
                rho=apply(rho,op;noise=nm)
            end
        # end
    end

    return rho
end

"""
`_final_measurement(state::AbstractVectorS, options::Options)`

Performs the final measurement on a quantum state based on specified options.

- `state`: The quantum state vector (sa.SparseVector) to be measured.
- `options`: Measurement options specifying the basis and error model.

Return the state vector to reflect the measurement outcome.
"""
function _final_measurement(state::AbstractVectorS,options::Options)#todo parallel measurement use apply_all()

    N=get_N(state)
    final_noise=options.readout_noise #measurement noise

    # final measurement_error or random measurement
    if options.measurement_basis=="Z"
        for qubit=1:N

            o=Op("ZBasis",gate.I,qubit)
            state=o.expand(N)*state

            if isa(final_noise,QuantumChannel)
                state=apply_noise(state, o, final_noise)
            end

        end
    elseif options.measurement_basis=="X"
        for qubit=1:N

            o=Op("XBasis",gate.H,qubit)
            state=o.expand(N)*state

            if isa(final_noise,QuantumChannel)
                state=apply_noise(state, o, final_noise)
            end

        end
    elseif options.measurement_basis=="Y"
        for qubit=1:N

            o=Op("YBasis",gate.HSp,qubit)
            state=o.expand(N)*state

            if isa(final_noise,QuantumChannel)
                state=apply_noise(state, o, final_noise)
            end

        end
    elseif options.measurement_basis=="R"
            # random measurement basis
        for qubit=1:N
            rOp=rand(1:3)

            o=Op("RBasis",[gate.H,gate.HSp,gate.I][rOp],qubit)
            state=o.expand(N)*state

            if isa(final_noise,QuantumChannel)
                state=apply_noise(state, o, final_noise)
            end

        end
    else
        throw("measurement_basis error!")
    end

    return state

end

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

