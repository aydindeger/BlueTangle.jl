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

"""
    AnsatzOptions(; N::Int, ops::Union{Vector{String},Vector{<:QuantumOps}}, noise=false, init::Union{sa.SparseVector,Circuit}=sa.sparse([]), model::String="lbfgs", number_of_iterations::Int=1000, learning_rate::Float64=0.01, pars_initial::Vector=[], deep_circuit::Bool=false, history::Bool=true)

Constructs an `AnsatzOptions` object that contains the configuration for a variational quantum circuit ansatz.

# Arguments
- `N::Int`: The number of qubits in the ansatz.
- `ops::Union{Vector{String},Vector{<:QuantumOps}}`: The quantum operations defining the ansatz, either as a vector of strings or a vector of `QuantumOps`.
- `noise=false`: Specifies whether to include noise in the ansatz. Can be either a `NoiseModel` or a boolean value.
- `init::Union{sa.SparseVector,Circuit}=sa.sparse([])`: The initial state of the ansatz, either as a sparse vector or a `Circuit`. Defaults to the zero state.
- `model::String="lbfgs"`: The optimization model to use. Can be "lbfgs", "adam", "descent", "radam", "momentum", or "nesterov".
- `number_of_iterations::Int=1000`: The number of iterations for the optimization.
- `learning_rate::Float64=0.01`: The learning rate for the optimization.
- `pars_initial::Vector=[]`: The initial parameters for the ansatz. If not provided, random parameters are generated.
- `deep_circuit::Bool=false`: Specifies whether to use a deep circuit ansatz.
- `history::Bool=true`: Specifies whether to record the optimization history.

# Returns
- An `AnsatzOptions` object containing the configuration for the variational quantum circuit ansatz.
"""
struct AnsatzOptions
    N::Int
    ops::Vector{<:QuantumOps}
    args::Vector{Int}
    loss::Function
    noise::Union{NoiseModel,Bool}
    dim::Int
    pars_initial::Vector
    init::Union{sa.SparseVector,it.MPS,Circuit}
    number_of_iterations::Int
    model::String
    learning_rate::Float64
    deep_circuit::Bool
    optimizer::Union{Optimisers.Leaf,OptimKit.LBFGS}
    history::Bool
    
    function AnsatzOptions(;
        N::Int,
        ops::Union{Vector{String},Vector{<:QuantumOps}},
        loss::Union{Function,sa.SparseMatrixCSC}, #input state returns number
        noise=false,
        init::Union{sa.SparseVector,it.MPS,Circuit}=sa.sparse([]),
        model::String="lbfgs",
        number_of_iterations::Int=1000,
        learning_rate::Float64=0.01,
        pars_initial::Vector=[],
        deep_circuit::Bool=false,
        history::Bool=true
        )

        if isa(ops,Vector{String})
            ops_final,op_args_list,dim=_variational_circuit_from_string(N,ops;deep_circuit=deep_circuit)
        else #find arg numbers
            if deep_circuit
                println("note that deep_circuit argument in invalid for given set of operations")
            end
            op_args_list=Vector{Int}()
            dim=0
            for o=ops

                if o.q!=1 && abs(o.qubit-o.target_qubit)>1
                    throw("non-local gate $(o.name) is not allowed! Use control parameter instead or add swaps")
                end
                
                if isa(o.mat,Function)
                    arg_no=BlueTangle._find_argument_number(o.mat)
                    dim += arg_no
                else
                    arg_no=0
                end
                push!(op_args_list,arg_no)
            end
            ops_final=ops
        end
    
        if init==sa.sparse([])
            state = zero_state(N)
        elseif typeof(init)==Circuit
            state=to_state(circuit)
        else
            state=init
        end

        if dim < length(pars_initial)
            println("Number of parameters do not match with the ansatz.\n$(dim) parameters will be used.")
            pars_initial=pars_initial[1:dim]
        else
            pars_initial=rand(dim) * pi
        end

        if lowercase(model)=="adam"
            optimizer=Optimisers.setup(Optimisers.Adam(learning_rate),pars_initial)
            # optimizer=Optimisers.setup(Optimisers.OptimiserChain(Adam(learning_rate)),pars_initial)
        elseif lowercase(model)=="descent" || lowercase(model)=="gradient"
            optimizer=Optimisers.setup(Optimisers.Descent(learning_rate),pars_initial)
        elseif lowercase(model)=="radam"
            optimizer=Optimisers.setup(Optimisers.RAdam(learning_rate),pars_initial)
        elseif lowercase(model)=="momentum"
            optimizer=Optimisers.setup(Optimisers.Momentum(learning_rate),pars_initial)
        elseif lowercase(model)=="nesterov"
            optimizer=Optimisers.setup(Optimisers.Nesterov(learning_rate),pars_initial)
        elseif lowercase(model)=="lbfgs" #OptimKit.jl LBFGS
            optimizer = OptimKit.LBFGS(; maxiter=number_of_iterations)
        else
            optimizer = OptimKit.LBFGS(; maxiter=number_of_iterations)
        end

        if isa(loss,sa.SparseMatrixCSC)
            loss_new(state)=isa(state,it.MPS) ? throw("loss is not MPO.") : real(state' * loss * state)
        else
            loss_new=loss
        end

        return new(N,ops_final,op_args_list,loss_new,noise,dim,pars_initial,state,number_of_iterations,model,learning_rate,deep_circuit,optimizer,history)
    end

end


function _variational_circuit_from_string(N::Int,ops::Vector{String};deep_circuit::Bool=false)
    
    ops_final=Vector{QuantumOps}()
    op_args_list=Vector{Int}()
    dim=0

    brick_c=0
    for gate_name=ops

        gate_name=BlueTangle._clean_name(gate_name)

        c=1
        mb=deep_circuit ? 0 : mod(brick_c,2)
        while c<=N #runs over qubits

            if BlueTangle._name_with_two_qubit_gates_bool(gate_name) #two qubits

                if mb+c>=N
                    break
                end

                if BlueTangle._name_with_arg_bool(gate_name)
                    
                    f=eval(Symbol("_"*gate_name))
                    arg_no=BlueTangle._find_argument_number(f)
                    
                    r=Op(gate_name,f,mb+c,mb+c+1)
                    dim += arg_no

                else
                    arg_no=0
                    r=Op(gate_name,mb+c,mb+c+1)
                end

                c=deep_circuit ? c+1 : c+2

            else #one qubit
                
                if BlueTangle._name_with_arg_bool(gate_name)
                    f=eval(Symbol("_"*gate_name))
                    arg_no=BlueTangle._find_argument_number(f)
                    r=Op(gate_name,f,c)
                    dim += arg_no
                else
                    arg_no=0
                    r=Op(gate_name,c)
                end

                c=c+1

            end

            push!(op_args_list,arg_no)
            push!(ops_final,r)

        end

        if BlueTangle._name_with_two_qubit_gates_bool(gate_name) #two qubits
            brick_c += 1
        end

    end

    return ops_final,op_args_list,dim

end


"""
    variational_apply(pars::Vector, opt::AnsatzOptions)

Applies the variational quantum circuit defined by the `AnsatzOptions` to the initial state using the given parameters.

# Arguments
- `pars::Vector`: The parameters for the variational circuit.
- `opt::AnsatzOptions`: The `AnsatzOptions` object containing the configuration for the variational circuit.

# Returns
- The final state after applying the variational circuit.
"""
function variational_apply(pars::Vector,opt::AnsatzOptions)

    state=opt.init
    N=isa(state,it.MPS) ? get_M(state) : opt.N
    noise=opt.noise

    if isa(noise, NoiseModel)

        c=1
        for (op,fn)=zip(opt.ops,opt.args)

            if op.q!=1 && abs(op.qubit-op.target_qubit)>1
                throw("non-local gate $(op.name) is not allowed! Use control parameter instead or add swaps")
            end

            state=fn>0 ? op.expand(N,pars[c:c+fn-1]...)*state : op.expand(N)*state
            c = c+fn

            selected_noise = op.q == 1 ? noise.q1 : noise.q2
            if isa(selected_noise, QuantumChannel) && op.noisy==true
                state = apply_noise(state, op, selected_noise)
            end
        end

    else

        c=1
        for (op,fn)=zip(opt.ops,opt.args)

            if op.q!=1 && abs(op.qubit-op.target_qubit)>1
                throw("non-local gate $(op.name) is not allowed! Use control parameter instead or add swaps")
            end

            state=fn>0 ? op.expand(N,pars[c:c+fn-1]...)*state : op.expand(N)*state
            c = c+fn
        end

    end

    return state
    
end


"""
    VQE(opt::AnsatzOptions)

Performs the Variational Quantum Eigensolver (VQE) algorithm to find the optimal state of a given loss function using the variational quantum circuit defined by the `AnsatzOptions`.

# Arguments
- `opt::AnsatzOptions`: The `AnsatzOptions` object containing the configuration for the variational quantum circuit.

# Returns
- A tuple containing:
    - The optimization history (energy values) if `history=true`, otherwise the final energy value.
    - The optimized parameters for the variational circuit.
    - The final state obtained from the optimized parameters.
"""
function VQE(opt::AnsatzOptions)

    function loss_func(pars::Vector, opt::AnsatzOptions)
        state = variational_apply(pars, opt)
        return opt.loss(state)
    end

    N=opt.N
    pars = opt.pars_initial
    # learning_rate=opt.learning_rate
    energy_history = Vector{Float64}()#(undef,en_size)
    optimizer = opt.optimizer

    if lowercase(opt.model)=="lbfgs"

        function loss_func_and_grad(pars::Vector, opt::AnsatzOptions)
            state = variational_apply(pars, opt)
            gr=ForwardDiff.gradient(p -> loss_func(p,opt), pars)
            return (opt.loss(state),gr)
        end

        pars, res, gs, niter, normgradhistory = OptimKit.optimize(p -> loss_func_and_grad(p, opt), pars, optimizer)

        if opt.history==false
            return res,pars,variational_apply(pars, opt)
        else
            return normgradhistory[:,1],pars,variational_apply(pars, opt)
        end

    else

        for i in 1:opt.number_of_iterations
            
            gr=ForwardDiff.gradient(p -> loss_func(p,opt), pars)

            # pars -= learning_rate * gr
            optimizer, pars = Optimisers.update(optimizer, pars, gr)

            if opt.history==true# && mod(i,Int(opt.number_of_iterations/100))==0
                push!(energy_history,loss_func(pars, opt))
            end

        end

        if opt.history==false
            return loss_func(pars, N),pars,variational_apply(pars, opt)
        else
            return energy_history,pars,variational_apply(pars, opt)
        end

    end

end

