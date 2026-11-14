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

"""
    nparameters(ops)

Count the consecutive parameters consumed by local function-valued gates.
Fixed gates, channels, and measurements consume no parameters.
Binding applies to the top-level operation list; keep trainable gates outside
named static blocks.
"""
nparameters(ops::AbstractVector{<:QuantumOps}) = sum(_parameter_count, ops; init=0)

const _prima_models = (:prima, :cobyla, :uobyqa, :newuoa, :bobyqa, :lincoa)
const _gradient_models = (:adam, :descent, :gradient, :radam, :momentum, :nesterov, :lbfgs)
_optimizer_name(model) = Symbol(lowercase(String(model)))

function _new_optimizer(model, parameters, iterations, learning_rate)
    model = _optimizer_name(model)
    model in _prima_models && return nothing
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
    AnsatzOptions(; N, ops, loss, model="cobyla", pars_initial=[], kwargs...)

Compatibility configuration for `VQA`/`VQE`. New workflows can use an operation
list, `variational_apply(parameters, ops)`, and `train(loss, parameters)` directly.
`ops` accepts explicit operations or gate-name layers. Initial states can be
vectors, density matrices, MPS states, or circuits. A matrix/MPO `loss` denotes
an expectation value; a function receives the final state.

`number_of_iterations=1000` sets the iteration limit (evaluation budget for
PRIMA), `learning_rate=0.01` controls first-order optimizers, and `history=true`
keeps the legacy history return. `paramshift=true` requires a matrix/MPO loss
and supported built-in gate functions. Use an explicit density matrix for
automatic differentiation through noise, and avoid differentiating sampled
measurement outcomes. MPS truncation is not differentiable by ForwardDiff;
use `train(...; gradient=:finite_difference)` for a deterministic truncated loss.

`noise=false`, `init=zero_state(N)`, `cutoff=-1`, and `deep_circuit=false` retain
the existing defaults. Longer initial parameter vectors are truncated for
compatibility; `nparameters(ops)` gives the exact required size for new code.
"""
struct AnsatzOptions{F,S,O}
    N::Int
    ops::Vector{<:QuantumOps}
    args::Vector{Int}
    loss::F
    noise::Union{NoiseModel,Bool}
    dim::Int
    pars_initial::Vector
    init::S
    number_of_iterations::Int
    model::String
    learning_rate::Float64
    cutoff::Number
    deep_circuit::Bool
    optimizer::Union{Optimisers.Leaf,OptimKit.LBFGS,Nothing}
    paramshift::Bool
    history::Bool
    observable::O
end

function AnsatzOptions(; N::Int, ops::Union{Vector{String},Vector{<:QuantumOps}},
                       loss::Union{Function,AbstractMatrix,it.MPO}, noise=false,
                       init=sa.sparse([]), model::String="cobyla",
                       number_of_iterations::Int=1000, learning_rate::Real=0.01,
                       cutoff::Real=-1, pars_initial::AbstractVector=[],
                       deep_circuit::Bool=false, paramshift::Bool=false,
                       history::Bool=true)
    N > 0 || throw(ArgumentError("N must be positive"))
    number_of_iterations > 0 || throw(ArgumentError("number_of_iterations must be positive"))
    isfinite(learning_rate) && learning_rate > 0 || throw(ArgumentError("learning_rate must be positive and finite"))
    name = _optimizer_name(model)
    name in (_prima_models..., _gradient_models...) || throw(ArgumentError("unknown optimizer: $model"))
    cutoff == -1 || (isfinite(cutoff) && cutoff >= 0) || throw(ArgumentError("cutoff must be -1 or nonnegative"))
    state = init isa AbstractVector && isempty(init) ? zero_state(N) : copy(init)
    state isa Union{AbstractVector,AbstractMatrix,it.MPS} || throw(ArgumentError("init must be a vector, density matrix, or MPS"))
    cutoff != -1 && !(state isa it.MPS) && throw(ArgumentError("cutoff applies only to MPS states"))
    loss isa it.MPO && !(state isa it.MPS) && throw(ArgumentError("an MPO loss requires an MPS initial state"))
    get_N(state) == N || throw(DimensionMismatch("initial state does not have N qubits"))
    loss isa AbstractMatrix && size(loss) != (2^N,2^N) && throw(DimensionMismatch("loss matrix must act on N qubits"))
    operations = ops isa Vector{String} ? _variational_circuit_from_string(N,ops; deep_circuit=deep_circuit)[1] : copy(ops)
    get_N(operations) <= N || throw(DimensionMismatch("an ansatz operation exceeds N qubits"))
    counts = Int[_parameter_count(op) for op in operations]
    dim = sum(counts)
    parameters = isempty(pars_initial) ? rand(dim)*pi : collect(float.(pars_initial))
    length(parameters) >= dim || throw(DimensionMismatch("the ansatz requires $dim initial parameters"))
    parameters = parameters[1:dim]
    all(isfinite, parameters) || throw(ArgumentError("initial parameters must be finite"))
    observable = loss isa Union{AbstractMatrix,it.MPO} ? loss : nothing
    paramshift && isnothing(observable) && throw(ArgumentError("parameter shifts require a matrix/MPO expectation loss; differentiate nonlinear losses with train instead"))
    paramshift && state isa it.MPS && cutoff > 0 && throw(ArgumentError("MPS truncation invalidates exact parameter shifts; use cutoff=0 or finite differences"))
    objective = isnothing(observable) ? loss : state -> real(expect(state, observable))
    optimizer = _new_optimizer(name, parameters, number_of_iterations, learning_rate)
    return AnsatzOptions(N, operations, counts, objective, noise, dim, parameters,
                         state, number_of_iterations, String(name), Float64(learning_rate),
                         cutoff, deep_circuit, optimizer, paramshift, history, observable)
end

"""
    EfficientSU2(N, reps=1, gates=["RY", "RZ"]; deep_circuit=true)

Return `reps+1` local rotation layers, separated by `reps` CX layers. Adjacent
CXs form a chain when `deep_circuit=true`; otherwise their bonds alternate in
disjoint layers. `nparameters(ops)` counts the resulting trainable angles.
"""
function EfficientSU2(N::Int, reps::Int=1, gates::Vector{String}=["RY", "RZ"];
                      deep_circuit::Bool=true)
    reps >= 0 || throw(ArgumentError("reps must be nonnegative"))
    names = String[]
    for layer in 0:reps
        append!(names, gates)
        layer < reps && push!(names, "CX")
    end
    return _variational_circuit_from_string(N, names; deep_circuit=deep_circuit)[1]
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
    return _variational_circuit_from_string(N, names; deep_circuit=deep_circuit)[1]
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
    counts = Int[_parameter_count(op) for op in operations]
    return operations, counts, sum(counts)
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
    state isa it.MPS && throw(ArgumentError("ForwardDiff cannot reliably differentiate MPS decompositions; use train(...; gradient=:finite_difference) or an untruncated expectation parameter shift"))
    return nothing
end

"""
    variational_apply(parameters, ops; init=nothing, N=get_N(ops), noise=false, kwargs...)
    variational_apply(parameters, options::AnsatzOptions; noise_override=nothing, kwargs...)

Apply consecutive parameters to an operation list and return the final state.
The default input is `zero_state(N)`; explicit vector, density-matrix, and MPS
inputs are supported and are not modified. `nparameters(ops)` gives the required
parameter count. Vector gates act directly on amplitudes without building
full-register operators. Noise and conditional branches use ordinary `apply`
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
    return _variational_apply(parameters, ops, Int[_parameter_count(op) for op in ops], state, noise; kwargs...)
end

function _variational_apply(parameters, ops, counts, state, noise; kwargs...)
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

function variational_apply(parameters::AbstractVector, opt::AnsatzOptions;
                           noise_override::Union{Nothing,Bool,NoiseModel}=nothing, kwargs...)
    noise = isnothing(noise_override) ? opt.noise : noise_override
    settings = opt.init isa it.MPS && opt.cutoff != -1 ? merge((cutoff=opt.cutoff,), (;kwargs...)) : (;kwargs...)
    opt.init isa it.MPS && opt.paramshift && opt.cutoff == -1 && (settings = merge((cutoff=0.0,), settings))
    return _variational_apply(parameters, opt.ops, opt.args, opt.init, noise; settings...)
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
Nesterov, and the legacy PRIMA methods. PRIMA uses an objective-evaluation budget.
Histories record accepted iterates for LBFGS, updated iterates for first-order
methods, and objective evaluations for PRIMA; `history=false` returns an empty
history. `loss` always corresponds to the returned parameters.
"""
function train(loss, initial_parameters::AbstractVector; optimizer=:lbfgs,
               iterations::Int=100, learning_rate::Real=0.01,
               gradient=:forward, history::Bool=true)
    result = _train(loss, initial_parameters; optimizer=optimizer, iterations=iterations,
                    learning_rate=learning_rate, gradient=gradient, history=history)
    return (loss=result.loss, parameters=result.parameters, history=result.history)
end

function _train(loss, initial_parameters; optimizer=:lbfgs, iterations=100,
                learning_rate=0.01, gradient=:forward, history=true,
                value_gradient=nothing, optimizer_state=nothing, parameter_history=false)
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
    isfinite(learning_rate) && learning_rate > 0 || throw(ArgumentError("learning_rate must be positive and finite"))
    all(p -> p isa Real && isfinite(p), initial_parameters) || throw(ArgumentError("initial parameters must be finite real numbers"))
    parameters = collect(float.(initial_parameters))
    model = _optimizer_name(optimizer)
    model in (_prima_models..., _gradient_models...) || throw(ArgumentError("unknown optimizer: $optimizer"))
    losses = Float64[]
    parameter_values = Vector{typeof(parameters)}()
    if isempty(parameters) || iterations == 0
        value = loss(parameters)
        history && push!(losses, value)
        parameter_history && push!(parameter_values, copy(parameters))
        return (loss=value, parameters=parameters, history=losses, parameter_history=parameter_values)
    end
    if model in _prima_models
        objective = p -> begin
            value = loss(p)
            history && push!(losses, value)
            value
        end
        parameters, info = getproperty(pr, model)(objective, parameters; maxfun=iterations)
        return (loss=info.fx, parameters=parameters, history=losses, parameter_history=parameter_values)
    end
    fg = isnothing(value_gradient) ? _value_gradient(loss, parameters, gradient) : value_gradient
    state = isnothing(optimizer_state) ? _new_optimizer(model, parameters, iterations, learning_rate) : optimizer_state
    if model == :lbfgs
        parameters, value, _, _, trace = OptimKit.optimize(fg, parameters, state)
        history && append!(losses, trace[:,1])
    else
        for _ in 1:iterations
            _, derivative = fg(parameters)
            state, parameters = Optimisers.update(state, parameters, derivative)
            history && push!(losses, loss(parameters))
            parameter_history && push!(parameter_values, copy(parameters))
        end
        value = history ? last(losses) : loss(parameters)
    end
    return (loss=value, parameters=parameters, history=losses, parameter_history=parameter_values)
end

"""
    VQA(options::AnsatzOptions)
    VQE(options::AnsatzOptions)

Compatibility entry points. With `history=false`, return `(loss, parameters)`.
With history enabled, LBFGS returns `(history, parameters)`, first-order methods
return `(history, parameters, parameter_history)`, and PRIMA returns
`(loss, parameters)`. New code can use `train` for a consistent named result.
"""
function VQA(opt::AnsatzOptions)
    result = _train(p -> _loss(p,opt), opt.pars_initial; optimizer=opt.model,
                    iterations=opt.number_of_iterations, learning_rate=opt.learning_rate,
                    history=opt.history, optimizer_state=opt.optimizer,
                    value_gradient=opt.paramshift ? p -> loss_and_grad_paramshift(p,opt) : nothing,
                    parameter_history=opt.history && _optimizer_name(opt.model) in _gradient_models && opt.model != "lbfgs")
    !opt.history || _optimizer_name(opt.model) in _prima_models ||
        return opt.model == "lbfgs" ? (result.history,result.parameters) : (result.history,result.parameters,result.parameter_history)
    return result.loss, result.parameters
end
"""`VQE(options)` is the energy-minimization alias of [`VQA`](@ref)."""
VQE(opt::AnsatzOptions) = VQA(opt)
_loss(parameters, opt) = opt.loss(variational_apply(parameters, opt))

function _parameter_shift_rule(op::Op, parameter::Int)
    name = _clean_name(op.name)
    gate_function = op.mat
    if gate_function isa _AdjointGate
        gate_function = gate_function.parent
        name = rstrip(name, '†')
    end
    get(_parameterized_gates, name, nothing) === gate_function ||
        throw(ArgumentError("Parameter shifts require a built-in gate function; use automatic differentiation for custom functions."))

    # A generator with two eigenvalues separated by frequency ω has the
    # two-point rule ω/2 [f(θ+π/(2ω)) - f(θ-π/(2ω))].
    if name == "SWAPA"
        return ((0.5, pi/2),)
    end
    three_eigenvalues = name in ("GIVENS", "RXY") || (name == "FSIM" && parameter == 1)
    frequency = name in ("RXX", "RYY") ? 2.0 : 1.0
    pauli_rotation = name in ("RX", "RY", "RZ", "RXX", "RYY", "RZZ") || (name == "U3" && parameter == 1)
    if three_eigenvalues || (op.control > 0 && pauli_rotation)
        # Eigenvalues {-r,0,r} produce frequencies r and 2r. These two
        # symmetric shifts differentiate both frequencies exactly.
        r = three_eigenvalues ? 1.0 : frequency/2
        return ((pi/(4r), r/2*(1 + inv(sqrt(2)))),
                (3pi/(4r), r/2*(inv(sqrt(2)) - 1)))
    end
    return ((pi/(2frequency), frequency/2),)
end

function loss_and_grad_paramshift(p::AbstractVector{<:Real}, opt::AnsatzOptions)
    # Parameter-shift differentiation assumes loss is an expectation value
    # (linear in the output density matrix), as in VQE.
    isnothing(opt.observable) && throw(ArgumentError("parameter shifts require an expectation loss"))
    ℓ = _loss(p, opt)
    g = zeros(eltype(p), length(p))
    base = copy(p)
    i = 1
    for (op, fn) in zip(opt.ops, opt.args)
        for parameter in 1:fn
            for (shift, weight) in _parameter_shift_rule(op, parameter)
                base[i] = p[i] + shift
                fp = _loss(base, opt)
                base[i] = p[i] - shift
                fm = _loss(base, opt)
                g[i] += weight * (fp - fm)
            end
            base[i] = p[i]
            i += 1
        end
    end
    return (ℓ, g)
end
