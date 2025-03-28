"""
`QuantumOps`

Abstract type representing quantum operations.
"""
abstract type QuantumOps end


function __calc_prob(state::AbstractVectorS,kraus_vec::Vector,qubit::Int)

    pA=partial_trace(state,qubit)
    
    return [real(la.tr(kraus * pA * kraus')) for kraus in kraus_vec]#probs

end

function __calc_prob(state::AbstractVectorS,kraus_vec::Vector,qubit::Int,target_qubit::Int)

    distance=abs(qubit-target_qubit)

    if distance==1
        pA=partial_trace(state,qubit,target_qubit)#local
    else
        pA=partial_trace(state,[qubit,target_qubit])#nonlocal
    end
    
    return [real(la.tr(kraus * pA * kraus')) for kraus in kraus_vec]#probs

end

function __QuantumChannel_new_apply(state::AbstractVectorS,kraus_ops::Vector,qubit::Int)
    N=get_N(state)
    ind=BlueTangle._weighted_sample(BlueTangle.__calc_prob(state,kraus_ops,qubit))
    return sa.normalize(hilbert(N,kraus_ops[ind],qubit)*state)
end

function __QuantumChannel_new_apply(state::AbstractVectorS,kraus_ops::Vector,qubit::Int,target_qubit::Int)
    N=get_N(state)
    ind=BlueTangle._weighted_sample(BlueTangle.__calc_prob(state,kraus_ops,qubit,target_qubit))
    return sa.normalize(hilbert(N,kraus_ops[ind],qubit,target_qubit)*state)
end


function __calc_prob3(state::AbstractVectorS,kraus_vec::Vector,first_qubit::Int)

    pA=partial_trace(state,[first_qubit,first_qubit+1,first_qubit+2])
    return [real(la.tr(kraus * pA * kraus')) for kraus in kraus_vec]#probs

end

function __QuantumChannel_new_apply3(state::AbstractVectorS,kraus_ops::Vector,first_qubit::Int)
    N=get_N(state)
    ind=BlueTangle._weighted_sample(BlueTangle.__calc_prob3(state,kraus_ops,first_qubit))
    return sa.normalize(hilbert3(N,kraus_ops[ind],first_qubit)*state)
end


function __QuantumChannel_new_apply(rho::sa.SparseMatrixCSC,kraus_ops::Vector,qubit::Int)
    new_rho=zero(rho)
    N=get_N(rho)
    for k=kraus_ops
        ek=hilbert(N,k,qubit)
        new_rho += ek * rho * ek'
    end
    return new_rho
end

function __QuantumChannel_new_apply(rho::sa.SparseMatrixCSC,kraus_ops::Vector,qubit::Int,target_qubit::Int)
    new_rho=zero(rho)
    N=get_N(rho)
    for k=kraus_ops
        ek=hilbert(N,k,qubit,target_qubit)
        new_rho += ek * rho * ek'
    end
    return new_rho
end

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
    name::String
    p::Float64
    kraus::Vector{AbstractMatrix}
    prob::Function
    apply::Function

    function QuantumChannel(model::String,kraus_ops::Vector{<:AbstractMatrix},p::Float64)

        if !is_valid_quantum_channel(kraus_ops)
            throw("not valid kraus operators: not CPTP!")
        end

        q=Int(log2(size(kraus_ops[1],1)))
        model=lowercase(model)

        if q==1

            new_prob(state::AbstractVectorS,first_qubit::Int)=__calc_prob(state,kraus_ops,first_qubit)
            new_apply(state::AbstractVectorS,first_qubit::Int)=__QuantumChannel_new_apply(state,kraus_ops,first_qubit)
            new_apply(rho::sa.SparseMatrixCSC,first_qubit::Int)=__QuantumChannel_new_apply(rho,kraus_ops,first_qubit)

            return new(q,model,p,kraus_ops,new_prob,new_apply)

        elseif q==2

            new_prob2(state::AbstractVectorS,first_qubit::Int,second_qubit::Int)=__calc_prob(state,kraus_ops,first_qubit,second_qubit)
            new_apply2(state::AbstractVectorS,first_qubit::Int,second_qubit::Int)=__QuantumChannel_new_apply(state,kraus_ops,first_qubit,second_qubit)
            new_apply2(rho::sa.SparseMatrixCSC,first_qubit::Int,second_qubit::Int)=__QuantumChannel_new_apply(rho,kraus_ops,first_qubit,second_qubit)

            return new(q,model,p,kraus_ops,new_prob2,new_apply2)

        else
            throw("Noise models are available only for 1 and 2 qubits!")
        end
    end

    function QuantumChannel(q::Int,model::String,p::Float64)
        model=lowercase(model)
        kraus_ops=q==1 ? BlueTangle.noise_model1(model,p) : BlueTangle.noise_model2(model,p)
        return QuantumChannel(model,kraus_ops,p)
    end

    function QuantumChannel(model::String,p::Float64)
        model=lowercase(model)
        kraus_ops=BlueTangle.noise_model1(model,p)
        return QuantumChannel(model,kraus_ops,p)
    end

end

const AbstractQuantumChannel = Union{QuantumChannel, Bool}

"""
    Create a NoiseModel
"""
struct NoiseModel
    q1::AbstractQuantumChannel
    q2::AbstractQuantumChannel

    function NoiseModel(n1::AbstractQuantumChannel, n2::AbstractQuantumChannel)
        if isa(n1, QuantumChannel) && isa(n2, QuantumChannel)
            if n1.q != 1 || n2.q != 2
                throw("size error in noise model (1|2)")
            end
        elseif isa(n1, QuantumChannel) && n1.q != 1
            throw("size error in noise model (1)")
        elseif isa(n2, QuantumChannel) && n2.q != 2
            throw("size error in noise model (2)")
        end
        return new(n1, n2)
    end

    NoiseModel(q1str::Vector,q2str::Vector)=NoiseModel(Noise1(q1str[1],q1str[2]),Noise2(q2str[1],q2str[2]))

    NoiseModel(q1str::String,val1::Float64,q2str::String,val2::Float64)=NoiseModel(Noise1(q1str,val1),Noise2(q2str,val2))

    NoiseModel(q1str::Vector,q2bool::Bool)=NoiseModel(Noise1(q1str[1],q1str[2]),false)
    NoiseModel(q1bool::Bool,q2str::Vector)=NoiseModel(false,Noise2(q2str[1],q2str[2]))

    NoiseModel(model::String,p::Float64)=NoiseModel(Noise1(model,p),Noise2(model,p))

end


"""
    OpQC=QuantumChannel but like Op

    OpQC(name::String,kraus_ops::Vector{<:AbstractMatrix},qubit::Int,target_qubit::Int=-1)
    OpQC(q::Int,model::String,p::Float64,qubit::Int,target_qubit::Int=-1)
"""
struct OpQC <: QuantumOps
    q::Int
    name::String
    qubit::Int
    target_qubit::Int
    control::Int
    kraus::Vector{AbstractMatrix}
    noisy::Bool
    type::String
    prob::Function
    apply::Function

    #the following is like quantum ops
    function OpQC(name::String,kraus_ops::Vector{<:AbstractMatrix},qubit::Int,target_qubit::Int=-1;type="")

        if !is_valid_quantum_channel(kraus_ops)
            throw("not valid kraus operators: not CPTP!")
        end

        q=Int(log2(size(kraus_ops[1],1)))
        name=lowercase(name)

        if q==1

            new_prob_op(state::AbstractVectorS)=__calc_prob(state,kraus_ops,qubit)
            new_apply_op(state::AbstractVectorS)=__QuantumChannel_new_apply(state,kraus_ops,qubit)
            new_apply_op(rho::sa.SparseMatrixCSC)=__QuantumChannel_new_apply(rho,kraus_ops,qubit)

            return new(q,name,qubit,target_qubit,-2,kraus_ops,false,type,new_prob_op,new_apply_op)

        elseif q==2

            new_prob2_op(state::AbstractVectorS)=__calc_prob(state,kraus_ops,qubit,target_qubit)
            new_apply2_op(state::AbstractVectorS)=__QuantumChannel_new_apply(state,kraus_ops,qubit,target_qubit)
            new_apply2_op(rho::sa.SparseMatrixCSC)=__QuantumChannel_new_apply(rho,kraus_ops,qubit,target_qubit)

            return new(q,name,qubit,target_qubit,-2,kraus_ops,false,type,new_prob2_op,new_apply2_op)

        elseif q==3

            # println("!!! note that 8x8 kraus operators will be applied to qubits: $(qubit),$(qubit+1),$(qubit+2) !!!")

            new_prob3_op(state::AbstractVectorS)=__calc_prob3(state,kraus_ops,qubit)
            new_apply3_op(state::AbstractVectorS)=__QuantumChannel_new_apply3(state,kraus_ops,qubit)

            return new(q,name,qubit,-1,-2,kraus_ops,false,type,new_prob3_op,new_apply3_op)

        else
            throw("Noise models are available only for up to 3 qubits!")
        end
    end

    function OpQC(model::String,p::Float64,qubit::Int,target_qubit::Int=-1)
        q=target_qubit>0 ? 2 : 1
        model=lowercase(model)
        kraus_ops=q==1 ? BlueTangle.noise_model1(model,p) : BlueTangle.noise_model2(model,p)
        return OpQC(model,kraus_ops,qubit,target_qubit;type="$(p)")
    end

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
custom_noise(q::Int,name_of_model::String,kraus::Vector{AbstractMatrix})=QuantumChannel(q,name_of_model,is_valid_quantum_channel(kraus) ? kraus : throw("define valid kraus operators"))

"""
`is_valid_quantum_channel(kraus::Vector{Matrix}) -> Bool`

Determines the validity of a set of Kraus operators.

- `kraus`: A vector of matrices, each representing a Kraus operator.

This function checks if the provided Kraus operators form a valid quantum channel. 
It does so by verifying if the sum of the products of each Kraus operator and its adjoint 
(approximately) equals the identity matrix. Returns `true` if the set is valid, `false` otherwise.
"""
function is_valid_quantum_channel(kraus_operators::Vector)
    n = size(kraus_operators[1],1)
    sum_kdagger_k = zeros(Complex{Float64}, n, n)
    
    # Iterate over each Kraus operator
    for K in kraus_operators
        sum_kdagger_k += K' * K
    end
    
    trace_preserving = sum_kdagger_k ≈ Matrix{Complex{Float64}}(la.I, n, n)
    
    choi_matrix = kraus_to_choi(kraus_operators)

    eigvals = la.eigen(choi_matrix).values

    completely_positive = (choi_matrix ≈ choi_matrix') && all(round.(eigvals,digits=10) .>= 0)#ishermitian and is_positive_semidefinite
    
    return trace_preserving && completely_positive
end

function kraus_to_choi(kraus_operators::Vector)
    d = size(kraus_operators[1], 1)
    choi = zeros(ComplexF64, d * d, d * d)
    for K in kraus_operators
        vK = vec(K)
        choi += vK * vK'
    end
    return choi
end
# function is_valid_quantum_channel(kraus::Vector)
#     sumk = sum(k' * k for k in kraus)
#     sumk ≈ Matrix(sa.I, size(sumk, 1), size(sumk, 1))
# end

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

##

"""
`Op(q::Int, name::String, mat::AbstractMatrix, qubit::Int, target_qubit::Int, noise::QuantumChannel) <: QuantumOps`

Represents a quantum operation.

- `q`: Number of qubits involved in the operation.
- `name`: Name of the operation.
- `mat`: Matrix representation of the quantum operation.
- `qubit`: index of the target qubit.
- `target_qubit`: index of the target qubit for two-qubit operations.
- `noise`: Noise model associated with the operation.

Constructs an Op object representing a quantum operation with optional noise.
"""
struct Op <: QuantumOps
    q::Int
    name::String
    qubit::Int
    target_qubit::Int
    control::Int
    mat::Union{AbstractMatrix, Function}
    noisy::Bool
    type::String
    expand::Function

    # Inner-constructor for gates defined from a function
    function Op(name::String,f::Function,qubit::Int,target_qubit::Int;type=nothing,noisy::Bool=true,control::Int=-2)
        if !isnothing(type)
            throw(ArgumentError("setting `type` for gates constructed from functions is not supported"))
        end

        q = _get_op_num_qubits(qubit, target_qubit, control)

        new_expand = _get_new_expand(f, qubit, target_qubit, control)

        return new(q,name,qubit,target_qubit,control,f,noisy,"f",new_expand)
    end
    # Inner-constructor for gates defined from a built-in or matrix 
    function Op(name::String,mat::AbstractMatrix,qubit::Int,target_qubit::Int;type::String="",noisy::Bool=true,control::Int=-2)

        q = _get_op_num_qubits(qubit, target_qubit, control)

        sizem = size(mat)
        if sizem != (2^q, 2^q)
            throw(ArgumentError("size of matrix $sizem not compatible with $(q)-qubit operation"))
        end

        if q == 1 && _is_it_measurement(name)

            if control!=-2
                throw(ArgumentError("measurement and control operations are incompatible."))
            end

            type = "🔬"
            noisy = false
            ismeasure = true
        else
            _name_with_arg_bool(name) ? type = "phase" : type="op"
            ismeasure = false
        end

        new_expand = _get_new_expand(name, mat, qubit, target_qubit, control, ismeasure)

        return new(q,name,qubit,target_qubit,control,mat,noisy,type,new_expand)
    end
end

# One qubit constructor with built-in gate
function Op(name::String,qubit::Int;kwargs...)

    if uppercase(name)=="RES"
        return ifOp("MZ",qubit,"I","X")
    else
        return Op(name, qubit, -1; kwargs...)
    end
end
# Two qubit constructor with built-in gate
function Op(name::String,qubit::Int,target_qubit::Int;kwargs...)

    mat=gates(name)

    return Op(name, mat, qubit, target_qubit; kwargs...)
end
# Three qubit constructor with built-in gate
function Op(name::String,qubit::Int,control_qubit::Int,target_qubit::Int)

    if name=="CCZ"
        return Op("CZ",qubit,target_qubit;control=control_qubit)
    elseif name=="CCX"
        return Op("CX",qubit,target_qubit;control=control_qubit)
    elseif name=="CCY"
        return Op("CY",qubit,target_qubit;control=control_qubit)
    elseif name=="CSWAP"
        return Op("SWAP",qubit,target_qubit;control=control_qubit)
    else
        throw("Unsupported three-qubit operation")
    end
end
# One qubit constructor with matrix or function
function Op(name::String,matf::Union{Function, AbstractMatrix},qubit::Int; kwargs...)
    return Op(name, matf, qubit, -1; kwargs...)
end

function _get_op_num_qubits(qubit::Int, target_qubit::Int, control::Int)
    if target_qubit == -1 
        # single qubit operator
        if qubit == control
            throw(ArgumentError("`qubit` must differ from `control` qubit"))
        else
            q = 1
        end
    else
        # two qubit operator
        if qubit == target_qubit
            throw(ArgumentError("`qubit` and `target_qubit` must differ"))
        elseif qubit == control && target_qubit == control
            throw(ArgumentError("either `qubit` or `target_qubit` must differ from `control` qubit"))
        else
            q = 2
        end
    end
    return q
end

function _get_new_expand(name::AbstractString, mat::AbstractMatrix, qubit::Int, target_qubit::Int, control::Int, ismeasure::Bool)
    function new_expand(N::Int)
        if target_qubit == -1
            if ismeasure
                rv = __measurement_hilbert(N,name,qubit)
            else
                rv = hilbert(N,mat,qubit;control=control)
            end
        else
            rv = hilbert(N,mat,qubit,target_qubit;control=control)
        end

        return rv
    end


    function new_expand(sites::Vector)
        if control != -2
            if uppercase(name) in ["X", "Y", "Z", "T", "TD", "S", "SD", "I", "H"]
                rv = _mat_to_tensor(sites, gate[Symbol("C$(uppercase(name))")], control, qubit)
            elseif name == "CX"
                rv = _mat_to_tensor(sites, gate.CCX, qubit, target_qubit; control=control)
            elseif name == "CZ"
                rv = _mat_to_tensor(sites, gate.CCZ, qubit, target_qubit; control=control)
            else
                throw("control tensor is not supported")
            end
            return rv
        end
    
        if target_qubit == -1
            rv = _mat_to_tensor(sites, mat, qubit)
        elseif target_qubit > 0
            size(mat, 1) == 4 ? (rv = _mat_to_tensor(sites, mat, qubit, target_qubit)) : throw("target_qubit is not defined")
        end
    
        return rv
    end

    return new_expand
end

function _get_new_expand(f::Function, qubit::Int, target_qubit::Int, control::Int)
    
    function new_expand(N::Int,pars...)
        if target_qubit == -1
            rv = hilbert(N,f(pars...),qubit;control=control)
        else
            rv = hilbert(N,f(pars...),qubit,target_qubit;control=control)
        end
        return rv
    end

    function new_expand(sites::Vector,pars...)
        # if (target_qubit > 0 && abs(qubit-target_qubit)>1) || control!=-2
        #     throw("nonlocal tensor is not supported")
        # end

        if control!=-2
            throw("nonlocal tensor is not supported")
        end

        if target_qubit == -1
            rv = _mat_to_tensor(sites,f(pars...),qubit)
        else
            rv = _mat_to_tensor(sites,f(pars...),qubit,target_qubit)
        end

        return rv
    end

    return new_expand
end

Op(nameAngle::Vector,qubit::Int;type::String="",noisy::Bool=true,control::Int=-2)=Op("$(nameAngle[1])($(nameAngle[2]...))",qubit;type=type,noisy=noisy,control=control)
Op(nameAngle::Vector,qubit::Int,target_qubit::Int;type::String="",noisy::Bool=true,control::Int=-2)=Op("$(nameAngle[1])($(nameAngle[2]))",qubit,target_qubit;type=type,noisy=noisy,control=control)


function _measurement_mat(name::String)
    uppercase_name=uppercase(name)
    if uppercase_name=="M(Z)" || uppercase_name=="MZ"
        return BlueTangle.gate.I
    elseif uppercase_name=="M(X)" || uppercase_name=="MX"
        return BlueTangle.gate.H
    elseif uppercase_name=="M(Y)" || uppercase_name=="MY"
        return BlueTangle.gate.HSP
    end
end

function __measurement_hilbert(N::Int,name::String,qubit::Int)
    if uppercase(name)=="MR" || uppercase(name)=="M(R)"
        rname=rand(["MX","MY","MZ"])
        # println("random mid-measurement applied in $(rname) basis")
        BlueTangle.hilbert(N,_measurement_mat(rname),qubit)
    else
        BlueTangle.hilbert(N,_measurement_mat(name),qubit)
    end
end

function __ifOp_apply(state::AbstractVectorS,name::String,qubit::Int,if0::Vector{Op},if1::Vector{Op},noise::Union{Bool,QuantumChannel})

    N=get_N(state)
    rotMat=__measurement_hilbert(N,name,qubit)
    state=rotMat*state #rotate
    state,ind=born_measure_Z(N,state,qubit)
    state=rotMat'*state #rotate back

    if_op_list=ind==0 ? if0 : if1
    for ifop=if_op_list
        state=BlueTangle.apply(state,ifop;noise=noise)
    end

    return state,ind

end

function __ifOp_apply(rho::sa.SparseMatrixCSC,name::String,qubit::Int,if0::Vector{Op},if1::Vector{Op},noise::Union{Bool,QuantumChannel})

    throw("error: fix this!")
    N=get_N(rho)

    if uppercase(name)=="MR" || uppercase(name)=="M(R)"
        println("Note that the density matrix will be calculated based on a single random basis.
        You might need to average over many density matrices.")
    end

    mat=__measurement_hilbert(N,name,qubit)
    rho = mat * rho * mat'; #rotate

    born_ops=[gate.P0, gate.P1]

    eb=BlueTangle.hilbert(N,born_ops[1],qubit)
    rho1 = eb * rho * eb'
    for ifop=if0
        rho1=BlueTangle.apply(rho1,ifop;noise=noise)
    end

    eb=BlueTangle.hilbert(N,born_ops[2],qubit)
    rho = eb * rho * eb'
    for ifop=if1
        rho=BlueTangle.apply(rho,ifop;noise=noise)
    end

    return rho+rho1

end

"""
`ifOp(q::Int, name::String, mat::AbstractMatrix, qubit::Int, if01::Tuple{Matrix,Matrix}, noise::QuantumChannel) <: QuantumOps`

Represents a conditional quantum operation used for mid-circuit born measurements. It is specifically designed for mid-circuit measurements in the X, Y, Z, or a random basis (R). Depending on the measurement outcomes (0 or 1), different gates specified in `if01` can be applied.

# Fields
- `q`: Integer representing the number of qubits involved in the operation.
- `name`: String indicating the name of the operation. Valid names are "MX", "MY", "MZ", or "MR", corresponding to operations in the X, Y, Z basis, or a random basis (R), respectively.
- `mat`: Matrix representing the quantum operation.
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
    mat::AbstractMatrix
    qubit::Int
    if01::Tuple{Vector{Op},Vector{Op}}
    type::String
    expand::Function #this stays for MR
    born_apply::Function

    function ifOp(name::String,qubit::Int,if0::Vector{Op},if1::Vector{Op})

        new_expand(N::Int)=__measurement_hilbert(N,name,qubit)
        new_born_apply(state::AbstractVectorS,noise=false)=__ifOp_apply(state,name,qubit,if0,if1,noise)
        new_born_apply(rho::sa.SparseMatrixCSC,noise=false)=__ifOp_apply(rho,name,qubit,if0,if1,noise)

        function new_expand(sites::Vector)
            throw("error: fix")
            # _mat_to_tensor(sites,_measurement_mat(name),qubit)
        end

        return new(1, name, gates(name), qubit, (if0,if1), _is_it_measurement(name) ?  "🔬" : throw("select MX or MY or MZ or MR basis."),new_expand,new_born_apply)
    end
    
    ifOp(name::String,qubit::Int,if0::Op,if1::Op) = ifOp(name,qubit,[if0],[if1])
    ifOp(name::String,qubit::Int,if0::Op,if1::Vector{Op}) = ifOp(name,qubit,[if0],if1)
    ifOp(name::String,qubit::Int,if0::Vector{Op},if1::Op) = ifOp(name,qubit,if0,[if1])

    ifOp(name::String,qubit::Int,if0::String="I",if1::String="I") = ifOp(name,qubit,[Op(if0,qubit)],[Op(if1,qubit)])
end

_find_argument_number(func::Function)=length(methods(func)[1].sig.parameters)-1

"""
f(state) -> return state
"""
struct OpF <: QuantumOps
    q::Int
    name::String
    apply::Function
    data::Union{Function, AbstractMatrixS, Vector{<:QuantumOps}}

    function OpF(name::String,f::Function)

        new_apply(state::Union{AbstractVectorS,it.MPS};kwargs...)=f(state;kwargs...)

        return new(1,name,new_apply,f)
    end

    function OpF(name::String,mat::AbstractMatrixS)

        new_apply_mat(state::AbstractVectorS)=mat*state

        return new(1,name,new_apply_mat,mat)
    end
    
    function OpF(name::String,ops::Vector{T}) where T <: QuantumOps

        function new_apply_mat2(state::Union{AbstractVectorS,it.MPS};kwargs...)

            if isa(state,it.MPS)
                for o=ops
                    state=apply(state,o;kwargs...)
                end
            else
                for o=ops
                    state=o*state
                end
            end

            return state
        end

        return new(1,name,new_apply_mat2,ops)
    end
    
    #fix: add noise
end

function _target_find(op::T) where T<:QuantumOps

    if isa(op,ifOp) || isa(op,OpF)
        return -1
    else
        return op.target_qubit
    end
end

function _control_find(op::T) where T<:QuantumOps
    if isa(op,ifOp) || isa(op,OpQC) || isa(op,OpF)
        return -2
    else
        return op.control
    end
end



##========== Struct ==========

"""
    Create Layout
"""
struct Layout
    N::Int
    layout::Union{Matrix, sa.SparseMatrixCSC}
    connectivity::Float64
    neighbors::Dict
    geometry::sa.SparseMatrixCSC
    noisy_swap::Bool
    fswap::Bool
    check::Function
    swap::Function

    function Layout(layout::Union{Matrix,Vector,sa.SparseMatrixCSC};noisy_swap::Bool=false,fswap::Bool=false)

        # layout=Int.(layout .> 0) #only 0 and 1

        if isa(layout,Matrix)
            layout=sa.sparse(layout)
        elseif isa(layout,Vector)
            layout=sa.sparse(hcat(layout...))
        end

        qubit_num_to_pos, pos_to_qubit_num = BlueTangle.enumerate_qubits(layout)

        N=length(qubit_num_to_pos)

        connectivity=BlueTangle._average_connectivity(layout)
        geometry=BlueTangle.illustrate_physical_qubits(layout, qubit_num_to_pos)

        neighbors_dict = BlueTangle.generate_neighbors_dict(qubit_num_to_pos, pos_to_qubit_num)

        new_check(qubit_num1::Int, qubit_num2::Int)=qubit_num1 ∈ neighbors_dict[qubit_num2]
        new_make_swaps(op::QuantumOps)=_make_swaps(op,neighbors_dict;noisy_swap=noisy_swap,fswap=fswap)

        println("Physical Qubits:")
        num_rows, num_cols = size(layout)
        for row in 1:num_rows
            for col in 1:num_cols
                s=geometry[row, col]
                print(s==0 ? "." : s, "\t")
            end
            println()
        end

        new(N,layout,connectivity,neighbors_dict,geometry,noisy_swap,fswap,new_check,new_make_swaps);
    end
end

##
##========== Circuit ==========

"""
`Options(circuit_name::String, measurement_basis::String, final_measurement_noise::QuantumChannel, Noise1::QuantumChannel, Noise2::QuantumChannel, twirl::Bool, noisy_swap::Bool, density_matrix::Bool)`

Represents configuration options for a quantum circuit.

- `circuit_name`: Name of the circuit.
- `measurement_basis`: Measurement basis used in the circuit.
- `final_measurement_noise`: Noise model for final measurement error.
- `noise`: Noise model for single-qubit and two-qubit operations.
- `twirl`: Boolean flag for twirling operations.
- `noisy_swap`: Boolean flag for swap errors.
- `density_matrix`: Boolean flag to indicate if density matrix should be calculated.

Constructs an Options object with specified settings for a quantum circuit.
"""
mutable struct Options
    circuit_name::String
    measurement_basis::String
    noise::Union{NoiseModel, Bool}
    twirl::Bool
    readout_noise::AbstractQuantumChannel
    measurement_mitigate::Bool
    density_matrix::Bool

    # Constructor with default values
    Options(;
        circuit_name::String="circuit", 
        measurement_basis::String="Z",
        noise=false,
        twirl::Bool=false,
        readout_noise=false,
        measurement_mitigate::Bool=false,
        density_matrix::Bool=false
    ) = new(circuit_name, uppercase(measurement_basis), noise, twirl, readout_noise, measurement_mitigate, density_matrix)
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
    stats_init::NamedTuple
    options::Options
    layout::Layout
    layers::Vector{Vector{<:QuantumOps}}
end



sa.sparse(circ) = sa.sparse(ComplexF64, circ)
function sa.sparse(::Type{T}, circ::Circuit) where {T<:Number}
    N = circ.stats.N
    return prod(o.expand(N) for o=reverse(vcat(circ.layers...)))
end

function sa.sparse(ops::Vector{QuantumOps})
    N=get_stats(ops).N
    return prod(o.expand(N) for o=reverse(ops))
end

hilbert(circ::Circuit)=sa.sparse(ComplexF64, circ) #to_state(circ)==hilbert(circ)*state

"""
`Measurement(bitstr::Union{Vector, UnitRange}, sample::Vector, expect::Vector, mag_moments::Vector, measurement_basis::String, number_of_experiment::Int, circuit_name::String, number_of_qubits::Int, density_matrix::sa.SparseMatrixCSC)`

Represents the result of quantum measurements.

- `bitstr`: Basis of measurement represented as integers or a range.
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
    bitstr::Union{Vector, UnitRange}
    sample::Vector #probability
    expect::Vector
    mag_moments::Vector
    measurement_basis::String
    number_of_experiment::Int
    circuit_name::String
    number_of_qubits::Int
    density_matrix::sa.SparseMatrixCSC
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
