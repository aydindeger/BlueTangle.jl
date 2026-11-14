"""
`hamiltonian(N::Int, string_of_ops::Vector, boundary::String="open")`

Constructs a Hamiltonian matrix for a quantum system with `N` qubits. The Hamiltonian is 
built based on the operators and their corresponding couplings specified in `string_of_ops`.

The `string_of_ops` should be an alternating array of coupling constants and operator 
strings. For example, `[.1, "Z,Z", .5, "X"]` implies a system with alternating couplings `.1` and `.5`, 
and operators "Z,Z" and "X".

The `boundary` parameter specifies the boundary conditions of the system. It can be either 
"open" or "periodic". In the case of "open" boundary conditions, interactions are only included 
for sites within the system size. For "periodic" boundary conditions, the system is treated as 
if it were wrapped around itself, allowing interactions that cross the end and start of the 
chain.

# Arguments
- `N::Int`: The number of qubits in the system.
- `string_of_ops::Vector`: An alternating vector of coupling constants and operator strings.
- `boundary::String`: The boundary condition of the system, either "open" or "periodic" (default is "open").

# Returns
- `SparseMatrix`: The Hamiltonian matrix representing the specified quantum system.

# Example
```julia
N = 6
string_of_ops = [-J, "Z,Z", -h, "X"]
H = hamiltonian(N, string_of_ops, "open")
```

This function iterates over the operator strings, applies each operator to the appropriate qubits 
based on the boundary conditions, and scales them by their corresponding coupling constants 
to construct the Hamiltonian.
"""
function hamiltonian(N::Int,string_of_ops::Vector,boundary::String="open")

N > 0 || throw(ArgumentError("N must be positive."))
_validate_hamiltonian_terms(string_of_ops,boundary)

couplings=[]
ops=[]
for (i,op)=enumerate(string_of_ops)
    if isodd(i)
        push!(couplings,op)
    else
        push!(ops,op)
    end
end

H=sa.spzeros(2^N,2^N)

for (id,op)=enumerate(ops)
    len_op=length(split(op,","))

    if boundary=="open"
        for site=1:N-(len_op-1)
            H += couplings[id]*expand_multi_op(op,[i for i=site:site+len_op-1],N)
        end
    elseif boundary=="periodic"
        for site=1:N
            H += couplings[id]*expand_multi_op(op,[mod1(i,N) for i=site:site+len_op-1],N)
        end
    end

end

return H

end

"""
    hamiltonian(rows_cols::Union{Tuple{Int64, Int64},Vector{Int64}}, string_of_ops::Vector, boundary::String="open")

Constructs a 2D Hamiltonian matrix for a quantum system with `rows * cols` qubits, potentially with double counting.

# Arguments
- `rows_cols`: A tuple or vector specifying the dimensions of the 2D lattice (rows, columns).
- `string_of_ops::Vector`: An alternating vector of coupling constants and operator strings.
- `boundary::String`: The boundary condition of the system, either "open" or "periodic" (default is "open").

# Returns
- `SparseMatrixCSC`: The Hamiltonian matrix representing the specified 2D quantum system.
"""
function hamiltonian(rows_cols::Union{Tuple{Int64, Int64},Vector{Int64}}, string_of_ops::Vector, boundary::String="open")
    length(rows_cols) == 2 && all(>(0), rows_cols) || throw(ArgumentError("Lattice dimensions must be two positive integers."))
    _validate_hamiltonian_terms(string_of_ops,boundary)
    rows, cols = rows_cols
    N = rows * cols
    couplings = string_of_ops[1:2:end]
    ops = string_of_ops[2:2:end]
    
    H = sa.spzeros(ComplexF64, 2^N, 2^N)
    
    for (id, op) in enumerate(ops)
        len_op = length(split(op, ","))
        
        # Horizontal connections
        for row in 1:rows
            if boundary == "open"
                for col in 1:cols-(len_op-1)
                    qubits = [(row-1)*cols + c for c in col:col+len_op-1]
                    H += couplings[id] * expand_multi_op(op, qubits, N)
                end
            elseif boundary == "periodic"
                for col in 1:cols
                    qubits = [(row-1)*cols + mod1(c, cols) for c in col:col+len_op-1]
                    H += couplings[id] * expand_multi_op(op, qubits, N)
                end
            end
        end
        
        # Vertical connections (if len_op > 1)
        if len_op > 1
            for col in 1:cols
                if boundary == "open"
                    for row in 1:rows-(len_op-1)
                        qubits = [((r-1)*cols + col) for r in row:row+len_op-1]
                        H += couplings[id] * expand_multi_op(op, qubits, N)
                    end
                elseif boundary == "periodic"
                    for row in 1:rows
                        qubits = [(mod1(r, rows)-1)*cols + col for r in row:row+len_op-1]
                        H += couplings[id] * expand_multi_op(op, qubits, N)
                    end
                end
            end
        end
    end
    
    return H
end

function _validate_hamiltonian_terms(terms::Vector, boundary::String)
    boundary in ("open", "periodic") || throw(ArgumentError("Boundary must be open or periodic."))
    iseven(length(terms)) || throw(ArgumentError("Hamiltonian terms must alternate coupling and operator string."))
    all(x -> x isa Number, terms[1:2:end]) || throw(ArgumentError("Hamiltonian couplings must be numbers."))
    all(x -> x isa String && !isempty(strip(x)), terms[2:2:end]) || throw(ArgumentError("Hamiltonian operators must be nonempty strings."))
    return nothing
end

_parameter_count(op::QuantumOps) = op isa Op && op.mat isa Function ? _find_argument_number(op.mat) : 0

const _prima_models = (:prima, :cobyla, :uobyqa, :newuoa, :bobyqa, :lincoa)
const _gradient_models = (:adam, :descent, :gradient, :radam, :momentum, :nesterov, :lbfgs)
_optimizer_name(model) = Symbol(lowercase(String(model)))

function _new_optimizer(model, parameters, iterations, learning_rate)
    model == :lbfgs && return OptimKit.LBFGS(; maxiter=iterations, verbosity=0)
    rule = if model == :adam
        Optimisers.Adam(learning_rate)
    elseif model in (:descent, :gradient)
        Optimisers.Descent(learning_rate)
    elseif model == :radam
        Optimisers.RAdam(learning_rate)
    elseif model == :momentum
        Optimisers.Momentum(learning_rate)
    elseif model == :nesterov
        Optimisers.Nesterov(learning_rate)
    else
        throw(ArgumentError("unknown optimizer: $model"))
    end
    return Optimisers.setup(rule, parameters)
end

"""
    EfficientSU2(N, reps=1, gates=["RY", "RZ"]; deep_circuit=true)

Return `reps+1` local rotation layers, separated by `reps` CX layers. Adjacent
CXs form a chain when `deep_circuit=true`; otherwise their bonds alternate in
disjoint layers.
"""
function EfficientSU2(N::Int, reps::Int=1, gates::Vector{String}=["RY", "RZ"];
                      deep_circuit::Bool=true)
    reps >= 0 || throw(ArgumentError("reps must be nonnegative"))
    names = String[]
    for layer in 0:reps
        append!(names, gates)
        layer < reps && push!(names, "CX")
    end
    return _variational_circuit_from_string(N, names; deep_circuit=deep_circuit)
end

"""
    generate_ansatz_circuit(N, reps=1, gates=["RY", "RZ", "CX"]; deep_circuit=true)

Repeat the supplied gate-name layers `reps` times. Bare rotation names make
trainable gates; explicit angles such as `"RY(0.3)"` remain fixed. Unlike
`EfficientSU2`, no final rotation layer is added automatically.
"""
function generate_ansatz_circuit(N::Int, reps::Int=1, gates::Vector{String}=["RY", "RZ", "CX"];
                                 deep_circuit::Bool=true)
    reps >= 0 || throw(ArgumentError("reps must be nonnegative"))
    names = repeat(gates, reps)
    return _variational_circuit_from_string(N, names; deep_circuit=deep_circuit)
end

function _variational_circuit_from_string(N::Int, names::Vector{String}; deep_circuit::Bool=false)
    N > 0 || throw(ArgumentError("N must be positive"))
    operations = QuantumOps[]
    two_qubit_layer = 0
    for raw_name in names
        name = uppercase(strip(raw_name))
        base = _clean_name(name)
        two_qubit = _name_with_two_qubit_gates_bool(name)
        trainable = _name_with_arg_bool(name) && !occursin('(', name)
        first_site = two_qubit && !deep_circuit ? 1 + mod(two_qubit_layer,2) : 1
        stride = two_qubit && !deep_circuit ? 2 : 1
        last_site = two_qubit ? N-1 : N
        for site in first_site:stride:last_site
            target = two_qubit ? site+1 : -1
            op = trainable ? Op(base, _parameterized_gates[base], site, target) : Op(name, site, target)
            push!(operations, op)
        end
        two_qubit_layer += two_qubit
    end
    return operations
end

_stochastic_operation(op::QuantumOps) = op isa Union{OpQC,ifOp} || (op isa Op && op.type == "🔬")
_stochastic_operation(op::OpF) = op.data isa Vector{<:QuantumOps} && any(_stochastic_operation,op.data)
_has_dual_numbers(x::AbstractArray) = eltype(x) <: Union{ForwardDiff.Dual,Complex{<:ForwardDiff.Dual}}
_has_dual_numbers(x) = false
_has_dual_gate(op::QuantumOps) = false
_has_dual_gate(op::Op) = op.mat isa AbstractMatrix && _has_dual_numbers(op.mat)
_has_dual_gate(op::OpF) = op.data isa Vector{<:QuantumOps} && any(_has_dual_gate, op.data)
_has_dual_gate(op::ifOp) = any(branch -> any(_has_dual_gate, branch), op.if01)

function _check_variational_gradient(state, differentiating, noise, stochastic)
    differentiating || return nothing
    if !(state isa AbstractMatrix) && (noise isa NoiseModel || stochastic)
        throw(ArgumentError("automatic differentiation through sampled noise/measurements is not an ensemble gradient; use a density-matrix initial state"))
    end
    state isa it.MPS && throw(ArgumentError("ForwardDiff cannot reliably differentiate MPS decompositions; use train(...; gradient=:finite_difference)"))
    return nothing
end

"""
    variational_apply(parameters, ops; init=nothing, N=get_N(ops), noise=false, kwargs...)

Apply consecutive parameters to an operation list and return the final state.
The default input is `zero_state(N)`; explicit vector, density-matrix, and MPS
inputs are supported and are not modified. Each function-valued local gate
consumes one parameter per argument; fixed operations consume none. Vector
gates act directly on amplitudes without building full-register operators.
Noise and conditional branches use ordinary `apply`
semantics. Differentiating sampled trajectories is rejected: use a density
matrix for an exact noisy ensemble, or define an explicit sampling estimator.
"""
function variational_apply(parameters::AbstractVector, ops::Vector{<:QuantumOps};
                           init=nothing, N::Int=isnothing(init) ? get_N(ops) : get_N(init),
                           noise::Union{Bool,NoiseModel}=false, kwargs...)
    N > 0 || throw(ArgumentError("supply N or an initial state for an empty circuit"))
    state = isnothing(init) ? zero_state(N) : init
    state isa Union{AbstractVector,AbstractMatrix,it.MPS} || throw(ArgumentError("invalid initial state"))
    get_N(state) == N || throw(DimensionMismatch("initial state does not have N qubits"))
    get_N(ops) <= N || throw(DimensionMismatch("an operation exceeds N qubits"))
    counts = Int[_parameter_count(op) for op in ops]
    dim = sum(counts)
    length(parameters) == dim || throw(DimensionMismatch("the circuit requires $dim parameters"))
    differentiating = _has_dual_numbers(parameters) || _has_dual_numbers(state) ||
                      any(_has_dual_gate, ops)
    _check_variational_gradient(state, differentiating, noise, any(_stochastic_operation, ops))
    state = state isa AbstractVector ? Vector(state) : copy(state)
    offset = 1
    for (op, count) in zip(ops, counts)
        _check_variational_gradient(state, _has_dual_numbers(state), noise, _stochastic_operation(op))
        parameters_for_gate = @view parameters[offset:offset+count-1]
        if state isa Vector && op isa Op && op.type != "🔬"
            isempty(kwargs) || _unsupported_keyword_check(kwargs, Symbol[])
            state = _variational_gate!(state, op, parameters_for_gate)
            # A gate function can capture Dual parameters from an outer loss.
            _check_variational_gradient(state, _has_dual_numbers(state), noise, false)
            noise isa NoiseModel && (state = apply_noise(state, op, noise))
        elseif count > 0
            state = apply(state, op; pars=parameters_for_gate, noise=noise, kwargs...)
        else
            state = apply(state, op; noise=noise, kwargs...)
        end
        offset += count
    end
    return state
end

# One derivative implementation is shared by all gradient optimizers. Reuse
# ForwardDiff's setup and obtain value + gradient in the same evaluation.
function _value_gradient(loss, parameters, gradient)
    if gradient == :forward
        result = ForwardDiff.DiffResults.GradientResult(parameters)
        config = ForwardDiff.GradientConfig(loss, parameters)
        return p -> begin
            ForwardDiff.gradient!(result, loss, p, config)
            # OptimKit retains old gradients, so return an independent vector.
            ForwardDiff.DiffResults.value(result), copy(ForwardDiff.DiffResults.gradient(result))
        end
    elseif gradient == :finite_difference
        return p -> begin
            g = similar(p)
            shifted = copy(p)
            for i in eachindex(p)
                step = cbrt(eps(eltype(p))) * max(one(eltype(p)), abs(p[i]))
                shifted[i] = p[i] + step
                plus = loss(shifted)
                shifted[i] = p[i] - step
                minus = loss(shifted)
                shifted[i] = p[i]
                g[i] = (plus-minus)/(2step)
            end
            loss(p), g
        end
    elseif gradient isa Function
        return p -> (loss(p), gradient(p))
    end
    throw(ArgumentError("gradient must be :forward, :finite_difference, or a gradient function"))
end

"""
    train(loss, initial_parameters; optimizer=:lbfgs, iterations=100,
          learning_rate=0.01, gradient=:forward, history=true)

Minimize a real scalar function of parameters. Returns the same named tuple
`(loss, parameters, history)` for every optimizer; inputs are not modified.
Use a closure with `variational_apply` for VQE, regression, or other quantum
learning tasks. No configuration object or dataset wrapper is required.
User-supplied loss functions and callbacks must be deterministic when using
automatic differentiation or finite differences.

`gradient=:forward` uses automatic differentiation; `:finite_difference` handles
deterministic, non-AD-compatible calculations such as truncated MPS evolution.
A gradient function is also accepted. Sampled noise requires a deliberately
chosen stochastic estimator; use density evolution for deterministic noisy
training. Supported optimizers are LBFGS, Adam, descent, RAdam, momentum,
Nesterov, and the derivative-free PRIMA methods. PRIMA uses an
objective-evaluation budget.
Histories record accepted iterates for LBFGS, updated iterates for first-order
methods, and objective evaluations for PRIMA; `history=false` returns an empty
history. `loss` always corresponds to the returned parameters.
"""
function train(loss, initial_parameters::AbstractVector; optimizer=:lbfgs,
               iterations::Int=100, learning_rate::Real=0.01,
               gradient=:forward, history::Bool=true)
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
    isfinite(learning_rate) && learning_rate > 0 || throw(ArgumentError("learning_rate must be positive and finite"))
    all(p -> p isa Real && isfinite(p), initial_parameters) || throw(ArgumentError("initial parameters must be finite real numbers"))
    parameters = collect(float.(initial_parameters))
    model = _optimizer_name(optimizer)
    model in (_prima_models..., _gradient_models...) || throw(ArgumentError("unknown optimizer: $optimizer"))
    losses = Float64[]
    if isempty(parameters) || iterations == 0
        value = loss(parameters)
        history && push!(losses, value)
        return (loss=value, parameters=parameters, history=losses)
    end
    if model in _prima_models
        objective = p -> begin
            value = loss(p)
            history && push!(losses, value)
            value
        end
        parameters, info = getproperty(pr, model)(objective, parameters; maxfun=iterations)
        return (loss=info.fx, parameters=parameters, history=losses)
    end
    fg = _value_gradient(loss, parameters, gradient)
    state = _new_optimizer(model, parameters, iterations, learning_rate)
    if model == :lbfgs
        parameters, value, _, _, trace = OptimKit.optimize(fg, parameters, state)
        history && append!(losses, trace[:,1])
    else
        for _ in 1:iterations
            _, derivative = fg(parameters)
            state, parameters = Optimisers.update(state, parameters, derivative)
            history && push!(losses, loss(parameters))
        end
        value = history ? last(losses) : loss(parameters)
    end
    return (loss=value, parameters=parameters, history=losses)
end
