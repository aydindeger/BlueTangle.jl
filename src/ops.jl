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

"""
    shadow(ops, number_of_experiment::Int; options=Options(), N=get_N(ops),
           layout=nothing, qubits=nothing, backend=:statevector, kwargs...)

Estimate a reduced density matrix using uniformly random local X, Y, Z
measurements. `qubits` selects the retained qubits (an integer or a collection
of distinct indices). The returned sparse matrix has size `2^k × 2^k` for
`k` selected qubits, in ascending qubit order, as in `partial_trace`.
Omitting `qubits` preserves the full-system estimate.

Each shot contributes the tensor product, over the selected qubits only, of
`(I + 3*(-1)^b*P)/2`, where `P` is the measured Pauli and `b` is its outcome.
Tracing out other qubits simply removes their trace-one factors; no full
density matrix is constructed when requesting a subsystem.

The default `backend=:statevector` uses a `2^N` statevector. Use `backend=:tensor`
(or `:mps`) to evolve and sample an MPS directly, without constructing a full
statevector. The MPS backend accepts `cutoff` and `maxdim`; truncation introduces
simulation error in addition to sampling error, and cost depends on the bond
dimensions. Both backends repeat preparation for every shot, including any
circuit noise or mid-circuit measurements.

With ideal readout, the estimator is unbiased for the simulated state, Hermitian
and trace one, but need not be positive at finite shot counts. Its eigenvalues are not probabilities
and cannot be used directly to estimate entropy. Readout noise and readout
mitigation are unsupported and raise an `ArgumentError`. The final
`measurement_basis` is replaced by independent random Pauli bases.
Operations run in their listed order; routing is performed only when `layout`
is supplied. `N` can include idle qubits, and omitted `qubits` selects all qubits.

```julia
ops = [Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
rho_13 = shadow(ops, 10_000; qubits=[1, 3]) # only a 4 × 4 estimate
rho_2 = shadow(ops, 10_000; qubits=2)       # only a 2 × 2 estimate
rho_13_mps = shadow(ops, 10_000; qubits=[1, 3], backend=:mps, cutoff=1e-12)
```
"""
function shadow(ops::Vector{<:QuantumOps}, number_of_experiment::Int;
                options::Options=Options(), N::Int=get_N(ops),
                layout::Union{Nothing,Layout}=nothing, qubits=nothing,
                backend::Union{Symbol,AbstractString}=:statevector, kwargs...)
    number_of_experiment > 0 || throw(ArgumentError("shots must be positive"))
    N > 0 || throw(ArgumentError("N must be positive; supply N explicitly for an empty operation list"))
    N >= get_N(ops) || throw(ArgumentError("N must include every operated qubit"))
    ops = isnothing(layout) ? ops : route(ops, layout)
    N = max(N, get_N(ops))
    backend_name = _resolve_backend(backend)
    backend_name in ("statevector", "tensor") ||
        throw(ArgumentError("shadow supports :statevector or :tensor (:mps) backends"))
    _unsupported_keyword_check(kwargs, backend_name == "tensor" ? [:cutoff, :maxdim] : Symbol[])
    keep = isnothing(qubits) ? collect(1:N) : qubits isa Integer ? [qubits] : collect(qubits)
    all(q -> q isa Integer && 1 <= q <= N, keep) && allunique(keep) ||
        throw(ArgumentError("qubits must contain distinct integer indices between 1 and $N"))
    keep = sort!(Int.(keep))
    options.readout_noise === false && !options.measurement_mitigate ||
        throw(ArgumentError("shadow requires ideal readout; readout_noise and measurement_mitigate are not supported"))

    rho_construct = sa.spzeros(ComplexF64, 2^length(keep), 2^length(keep))
    rotations = (gate.H, gate.HSP, gate.I) # U† Z U = X, Y, Z; HSP = H*S†.
    local_estimators = [sa.sparse((gate.I + 3*(1-2*b)*P)/2)
                        for P in (gate.X, gate.Y, gate.Z), b in 0:1]
    measurement_basis = Vector{Int}(undef, N)
    for _ in 1:number_of_experiment
        # Fresh trajectory in the selected representation, including preparation noise.
        initial = backend_name == "tensor" ? zero_state(N_MPS(N)) : zero_state(N)
        state = apply(initial, _execution_ops(ops, options); noise=options.noise, kwargs...)
        for qubit in 1:N
            basis = rand(1:3)
            state = apply(state, Op(("Xbasis", "Ybasis", "Zbasis")[basis], rotations[basis], qubit); kwargs...)
            measurement_basis[qubit] = basis
        end
        # MPS sampling contracts tensors directly; joint bits preserve correlations
        # and also avoid packing long registers into a machine integer.
        bits = only(sample_bit(state, 1))

        # Average products within each shot, not products of single-qubit averages.
        snapshot = sa.sparse(ones(ComplexF64, 1, 1))
        for q in keep
            snapshot = kron(snapshot, local_estimators[measurement_basis[q], bits[q]+1])
        end
        rho_construct += snapshot
    end
    return rho_construct / number_of_experiment
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
    get_layers(ops::Vector{<:QuantumOps})

Return parallel layers, preserving operation order on every shared qubit.
Conditional operations reserve all branch qubits, and `OpF` separates layers.
Use this to inspect or select layers; `plotq(ops)` computes them automatically.
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

function _readout_response_inverse(readout_noise::QuantumChannel)
    readout_noise.q == 1 || throw(ArgumentError("readout noise must be a single-qubit channel"))
    response = zeros(Float64, 2, 2)
    for (outcome, projector) in enumerate((gate.P0, gate.P1))
        # Tr(E_m * rho) depends only on populations exactly when E_m is diagonal.
        # Sum first: coherence contributions from different Kraus terms can cancel.
        effect = sum(K' * projector * K for K in readout_noise.kraus)
        isapprox(effect, la.Diagonal(la.diag(effect)); atol=1e-12, rtol=1e-12) ||
            throw(ArgumentError("readout mitigation requires measured populations to be independent of input coherences"))
        response[outcome, :] = real.(la.diag(effect))
    end
    la.rank(response) == 2 ||
        throw(ArgumentError("readout mitigation requires an invertible classical response matrix"))
    return inv(response)
end

# Independent readout channels have a tensor-product classical response matrix.
function _measurement_mitigate_inv_matrix(N::Int, readout_noise::QuantumChannel)
    inverse = _readout_response_inverse(readout_noise)
    return foldl(kron, (inverse for _ in 1:N); init=ones(1, 1))
end

function _mitigate_probabilities(bitstr, probabilities, N, readout_noise)
    inverse = _readout_response_inverse(readout_noise)
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
    measure(state; label="state to measurement")
    measure(state, shots; label="state to measurement")

Return measurement statistics for a statevector, density matrix, or MPS.
Omitting `shots` uses exact probabilities; a positive count
uses sampled frequencies. Exact MPS enumeration is limited to 16 qubits, so
pass a shot count for larger states. Use `sample(state, shots)` for raw outcomes.
Large sampled MPS results carry an empty density-matrix placeholder.
"""
function measure(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, number_of_experiment::Int=-1; label="state to measurement")
    number_of_experiment == -1 || number_of_experiment > 0 || throw(ArgumentError("shots must be positive, or -1 for exact probabilities"))
    N = get_N(state)
    bitstr, probabilities = number_of_experiment == -1 ? sample_exact(state) : sample_state(state, number_of_experiment)
    rho = state isa AbstractMatrix ? sa.sparse(state) : _measurement_rho_placeholder(N)
    return _measurement_from_probs(bitstr, probabilities, N, number_of_experiment, "0", label, rho)
end

"""
    Measurement(samples::Vector{<:Integer}, N::Int; label="sample to measurement")

Summarize existing integer outcomes from an `N`-qubit register. Outcomes must
lie in `0:2^N-1`. This replaces `measure(samples, N)`, whose arguments could be
confused with an integer-valued statevector and a shot count.
"""
function Measurement(samples::Vector{<:Integer}, N::Int; label="sample to measurement")
    N >= 0 || throw(ArgumentError("N must be nonnegative"))
    isempty(samples) && throw(ArgumentError("at least one sample is required"))
    limit = big(1) << N
    all(outcome -> 0 <= outcome < limit, samples) ||
        throw(ArgumentError("samples must be integer outcomes in 0:2^N-1"))
    bitstr, probabilities = get_probs_from_sample(samples, N)
    return _measurement_from_probs(bitstr, probabilities, N, length(samples), "0", label, _measurement_rho_placeholder(N))
end

_is_cnot_like(op::QuantumOps) = op isa Op && op.q == 2 && op.control == -2 &&
    uppercase(op.name) in ("CNOT", "CX") && op.mat isa AbstractMatrix && op.mat == gate.CX

function _samples_to_measurement(samples::Vector{<:Integer}, N::Int, options::Options;
                                 rho=_measurement_rho_placeholder(N))
    bitstr, probabilities = get_probs_from_sample(samples, N)
    if options.measurement_mitigate && options.readout_noise isa QuantumChannel
        bitstr, probabilities = _mitigate_probabilities(bitstr, probabilities, N, options.readout_noise)
    end
    return _measurement_from_probs(bitstr, probabilities, N, length(samples), options.measurement_basis, options.circuit_name, rho)
end

function _execution_ops(ops::Vector{<:QuantumOps}, options::Options)
    options.twirl || return ops
    operations = QuantumOps[]
    for op in ops
        if _is_cnot_like(op)
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
    measure(ops, shots; options=Options(), N=get_N(ops), layout=nothing,
            backend=:statevector, kwargs...)

Return final-outcome statistics from fresh circuit shots. `sample(ops, shots)`
returns the corresponding raw integer outcomes. Choose statevector, MPS
(`:tensor`), or density-matrix (`:density_matrix`) evolution; MPS accepts `cutoff`
and `maxdim`. `N` can include idle qubits. Operations run in their listed order,
with routing only when `layout` is supplied. Readout mitigation and optional
density-matrix storage are configured through `Options`.
"""
function measure(ops::Vector{<:QuantumOps}, number_of_experiment::Int;
                 options::Options=Options(), N::Int=get_N(ops),
                 layout::Union{Nothing,Layout}=nothing,
                 backend::Union{Symbol,AbstractString}=:statevector, kwargs...)
    number_of_experiment > 0 || throw(ArgumentError("shots must be positive"))
    result = _sample_circuit(ops, number_of_experiment; options, N, layout, backend, kwargs...)
    return _samples_to_measurement(result.samples, result.N, options; rho=result.rho)
end

"""
    sample(ops, shots; options=Options(), N=get_N(ops), layout=nothing,
           backend=:statevector, kwargs...)

Return raw final outcomes as integers, starting each shot in `|0…0⟩`.
Qubit 1 is the most significant bit. Noisy statevector and MPS circuits use
independent trajectories; density-matrix evolution sums unobserved branches.
The circuit's explicit measurement records are available through `run`.

The backend, routing, measurement basis, readout noise, and twirling options
match `measure`. Use `measure` for summary statistics, readout mitigation, or
stored density matrices: raw outcomes cannot include those additional results.
"""
function sample(ops::Vector{<:QuantumOps}, shots::Int;
                options::Options=Options(), N::Int=get_N(ops),
                layout::Union{Nothing,Layout}=nothing,
                backend::Union{Symbol,AbstractString}=:statevector, kwargs...)
    options.measurement_mitigate &&
        throw(ArgumentError("raw samples cannot contain readout-mitigated weights; use measure(ops, shots; options=...)"))
    options.density_matrix &&
        throw(ArgumentError("raw samples cannot store a density matrix; use measure(ops, shots; options=...)"))
    return _sample_circuit(ops, shots; options, N, layout, backend, kwargs...).samples
end

function _sample_circuit(ops::Vector{<:QuantumOps}, number_of_experiment::Int;
                         options::Options, N::Int, layout::Union{Nothing,Layout},
                         backend::Union{Symbol,AbstractString}, kwargs...)
    number_of_experiment >= 0 || throw(ArgumentError("shots must be nonnegative"))
    N > 0 || throw(ArgumentError("N must be positive; supply N explicitly for an empty operation list"))
    N >= get_N(ops) || throw(ArgumentError("N must include every operated qubit"))
    ops = isnothing(layout) ? ops : route(ops, layout)
    N = max(N, get_N(ops))
    backend_name = _resolve_backend(backend)
    _unsupported_keyword_check(kwargs, backend_name == "tensor" ? [:cutoff, :maxdim] : Symbol[])
    T = N < 8sizeof(Int) ? Int : BigInt
    all_sample = T[]
    rho = _measurement_rho_placeholder(N)
    if options.density_matrix
        N <= 16 || throw(ArgumentError("density_matrix=true explicitly materializes 4^N entries; disable it for large MPS circuits"))
    end

    stochastic = options.noise != false || options.readout_noise != false ||
                 any(op -> op isa OpQC || op isa ifOp || op isa OpF || op.type == "🔬", ops)
    # Density evolution already sums noise and measurement branches exactly.
    reuse_state = (backend_name == "density_matrix" || !stochastic) &&
                  !any(_has_state_callback, ops) && !options.twirl && options.measurement_basis != "R"
    repetitions = reuse_state ? min(1, number_of_experiment) : number_of_experiment
    for _ in 1:repetitions
        initial = backend_name == "tensor" ? zero_state(N_MPS(N)) : zero_state(N)
        backend_name == "density_matrix" && (initial = initial * initial')
        state = apply(initial, _execution_ops(ops, options); noise=options.noise, kwargs...)
        state = _final_measurement(state, options; kwargs...)
        append!(all_sample, sample(state, reuse_state ? number_of_experiment : 1))
        if options.density_matrix
            vector = state isa it.MPS ? to_state(state) : state
            rho += (vector isa AbstractMatrix ? vector : vector * vector') / repetitions
        end
    end
    return (samples=all_sample, N=N, rho=sa.sparse(rho))
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

"""
    run(state, ops; noise=false, kwargs...)

Execute operations once and return `(state=final_state, record=bits)`.
The input state selects the backend. Operations run in their listed order;
`record` contains each explicit measurement outcome, including measurements
inside the selected feedback branch. Outcomes 0 and 1 denote the +1 and −1
Pauli eigenvalues. Reset adds no bit, and no final measurement is appended.

Density matrices are conditioned on the recorded measurement outcomes, while
noise and reset still use their complete Kraus sums. Statevectors and MPS states
follow one stochastic trajectory. MPS accepts `cutoff` and `maxdim`.
Use `apply(state, ops)` when only the evolved state is needed; density-matrix
`apply` averages over unobserved measurement outcomes.
"""
function run(state::Union{AbstractVector{<:Number},AbstractMatrixS,it.MPS},
             ops::Union{QuantumOps,Vector{<:QuantumOps},AbstractVector{<:Tuple}};
             noise::Union{NoiseModel,Bool}=false, kwargs...)
    haskey(kwargs, :track_measurements) &&
        throw(ArgumentError("run always returns the measurement record; omit track_measurements"))
    final_state, bits = apply(state, ops; noise, track_measurements=true, kwargs...)
    return (state=final_state, record=bits)
end

"""
    run(ops, shots=1; options=Options(), N=get_N(ops), layout=nothing,
        backend=:statevector, kwargs...)
    run(qasm_or_path::AbstractString, shots::Int=1; backend=:statevector, kwargs...)

Return a batch of measurement records, starting each shot in `|0…0⟩`.
This existing convenience form returns only records; `run(initial, ops)` returns
both the final state and record for one execution. Records follow operation-list
order and include measurements inside the chosen conditional branch. Reset adds
no bit, and no final measurement is appended; use `sample` for raw final outcomes
or `measure` for their statistics.

Choose `:statevector`, `:density_matrix`, or `:tensor` (`:mps`). MPS accepts
`cutoff` and `maxdim`, including for conditional branch operations. Density-matrix
execution conditions on the recorded measurements but still sums noise and
reset Kraus branches. `N` can include idle qubits. Routing is performed only when
`layout` is supplied. Final-readout options do not affect these explicit records.
"""
function run(ops::Vector{<:QuantumOps},number_of_experiment::Int=1;
             options::Options=Options(), N::Int=get_N(ops),
             layout::Union{Nothing,Layout}=nothing,
             backend::Union{Symbol,AbstractString}=:statevector,kwargs...)
    number_of_experiment > 0 || throw(ArgumentError("`number_of_experiment` must be positive."))
    N > 0 || throw(ArgumentError("N must be positive; supply N explicitly for an empty operation list"))
    N >= get_N(ops) || throw(ArgumentError("N must include every operated qubit"))
    ops = isnothing(layout) ? ops : route(ops, layout)
    N = max(N, get_N(ops))
    backend_name = _resolve_backend(backend)
    _unsupported_keyword_check(kwargs, backend_name == "tensor" ? [:cutoff, :maxdim] : Symbol[])
    measurements = Vector{Vector{Int}}(undef,number_of_experiment)
    for shot=1:number_of_experiment
        initial = backend_name == "tensor" ? zero_state(N_MPS(N)) : zero_state(N)
        backend_name == "density_matrix" && (initial = initial * initial')
        measurements[shot] = run(initial, _execution_ops(ops, options);
                                 noise=options.noise, kwargs...).record
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

function run(qasm_or_path::AbstractString,number_of_experiment::Int=1;
             options::Options=Options(),layout::Union{Nothing,Layout}=nothing,
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

For example, `get_operator("X", 3)` equals `kron(gate.X, gate.X, gate.X)`.
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
    string_to_matrix(list_of_operators::String)

Convert a string of comma-separated operators into a matrix representation.

- `list_of_operators`: Comma-separated gate names, such as `"Z,Z,RX(0.2)"`.

Returns a matrix representation of the operators.

```julia
matrix = string_to_matrix("Z,Z,RX(0.2)")
@assert matrix ≈ kron(gate.Z, gate.Z, _RX(0.2))
```
"""
function string_to_matrix(list_of_operators::String)
    ops_str=String.(split(list_of_operators,","));
    return foldl(kron,sa.sparse.(gates.(ops_str)))
end

"""
    hamming_distance(v1::Vector{Int}, v2::Vector{Int})
"""
hamming_distance(v1::Vector{Int}, v2::Vector{Int})=sum(v1 .⊻ v2)
