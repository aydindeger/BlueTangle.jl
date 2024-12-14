# state0=[1;0]
# state1=[0;1]

"""
`gate`

A constant that holds common quantum gates and projectors.

This includes:
- Identity (`I`)
- Pauli gates (`X`, `Y`, `Z`)
- Hadamard (`H`)
- Phase gates (`S`, `T`)
- Special gates like `sqrt(X)` (equal to RX(pi/2)*exp(1im*pi/4))
- Projectors (`P0` for |0><0|, `P1` for |1><1|)
- Controlled gates such as `CX` (CNOT), `CNOT` (an alias for CX), and `CZ`
- The `SWAP`, `ISWAP`, `FSWAP` gate
- The `ECR`, `SYC` gate

Each single-qubit gate is represented as a 2x2 matrix, while multi-qubit gates like `CNOT`, `ECR`, `SYC`, `CZ`, and `SWAP`, `ISWAP`, `FSWAP` are represented as 4x4 matrices.
"""
const gate = (
    I = [1. 0; 0 1],
    X  = [0 1.; 1 0],
    SX  = (1/2)*[1+1im 1.0-1im; 1-1im 1+1im],# SX ==RX(pi/2)*exp(1im*pi/4)
    XSQRT  = (1/2)*[1+1im 1.0-1im; 1-1im 1+1im],# SX ==RX(pi/2)*exp(1im*pi/4)
    Y  = [0 -im; im 0],
    Z  = [1.0 0; 0 -1],
    H  = (1/sqrt(2)) * [1.0 1; 1 -1],
    S  = [1. 0; 0 im],
    SD  = [1. 0; 0 -im], #S†
    T  = round.([1. 0; 0 exp(im * π / 4)],sigdigits=10),
    TD  = round.([1. 0; 0 exp(-im * π / 4)],sigdigits=10),
    HSP = (1/sqrt(2)) * [1.0+0im  0-1im;1+0im  0+1im],
    HY = (1/sqrt(2)) * [1+0.0im  0.0+1im; 1+0.0im  0.0-1im], #H*S
    H2 = [0.5 0.5 0.5 0.5;0.5 -0.5 0.5 -0.5;0.5 0.5 -0.5 -0.5;0.5 -0.5 -0.5 0.5],
    # HS = (1/sqrt(2)) * [1.0  1;1im  -1im],
    P0=[1.0 0; 0 0], #|0><0> #projector
    P1=[0 0; 0 1.0], #|1><1> #projector #n
    SP=[0 1;0 0],
    SM=[0 0;1 0],
    CX = [1. 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], # CNOT
    CNOT = [1. 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], # CNOT
    CY = [1. 0  0   0; 0 1.  0   0; 0 0  0 -im; 0 0 im  0],
    CZ = [1. 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1],
    SWAP =[1. 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1],
    ISWAP =[1. 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1],
    FSWAP =[1. 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1],
    SYC =[1. 0 0 0; 0 0 -im 0; 0 -im 0 0; 0 0 0 exp(-im*pi/6)],#Sycamore
    ECR = (1/sqrt(2)) * [0 + 0im 1. + 0im 0 + 0im 0 + 1im; 1 + 0im 0 + 0im 0 - 1im 0 + 0im; 0 + 0im 0 + 1im 0 + 0im 1 + 0im; 0 - 1im 0 + 0im 1 + 0im 0 + 0im], # (IX-XY)/sqrt(2),
    CCX = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0],
    CCZ = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0]
    )


one_qubit_gates=["I","X","Y","Z","SX","XSQRT","H","T","S","SD","P","U2","U3"]
two_qubit_gates=["CX","CNOT","CY","CZ","CP","RXX","GIVENS","FSIM","SWAP","ISWAP","FSWAP","SYC","ECR"]
gates_with_phase=["P","RX","RY","RZ","U2","U3","CP","GIVENS","FSIM","SWAPA","RXX","RYY","RZZ","RXY"]

    
"""
`random_gate_1(N::Int) -> Op`

Generates a random single-qubit gate operation.

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated single-qubit gate.
"""
function random_gate_1(N::Int)## N is total qubit number
    
    lambda=round(rand(),digits=2)
    op_name=rand(["I","X","Y","Z","H","S","SX","T","P","RX","RY","RZ"])

    if op_name in ["P","RX","RY","RZ"]
        name_return="$(op_name)($(lambda)π)"
    else
        name_return=op_name
    end

    return Op(name_return,gates(name_return),rand(1:N))
end

"""
`random_gate_2(N::Int) -> Op`

Generates a random two-qubit gate operation from ["CX","CZ","CP","GIVENS","FSIM","SWAP","ISWAP","FSWAP","SYC","ECR"]

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated two-qubit gate.
"""
function random_gate_2(N::Int)## N is total qubit number
    
    p1=round(rand(),digits=2)
    p2=round(rand(),digits=2)
    op_name=rand(two_qubit_gates)

    if op_name ∈ ["CP", "GIVENS","SWAPA","RXX","RYY","RZZ","RXY"]
        name_return="$(op_name)($(p1)π)"
    elseif op_name=="FSIM"
        name_return="$(op_name)($(p1)π,$(p2)π)"
    else
        name_return=op_name
    end

    r=rand(1:N-1)

    return Op(name_return,gates(name_return),r,r+1)
end

"""
`random_gate(N::Int) -> Op`

Generates a random quantum gate operation (either single or two-qubit).

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated quantum gate.
"""
random_gate(N)=rand([random_gate_1(N),random_gate_2(N)])

"""
`random_ops_2(N::Int, len::Int; measure_prob::Float64=0.0, measure_basis::Vector{String}=["MX", "MY", "MZ"]) -> Vector{QuantumOps}`

Create a sequence of random quantum gate operations, with optional mid-circuit measurements.

# Arguments
- `N::Int`: The number of qubits in the system.
- `len::Int`: The length of the sequence of operations to generate.

# Keyword Arguments
- `measure_prob::Float64`: The probability of adding a measurement operation after each gate.
- `measure_basis::Vector{String}`: The basis in which measurements are performed.

# Returns
- `Vector{QuantumOps}`: A vector of randomly chosen quantum operations (`QuantumOps`), each representing a gate or a measurement operation.

# Description
This function creates a vector of quantum operations, where each operation is either a randomly chosen gate from the set {"X", "Y", "Z", "H", "S", "CX","CZ","CP","GIVENS","FSIM","SWAP","ISWAP","FSWAP","SYC","ECR"} or a measurement operation, based on `measure_prob`. 

# Example
```julia
ops = random_ops(5, 10; measure_prob=0.2, measure_basis=["MX","MZ"])
```
This example generates a sequence of 10 random gates and measurements (with a 20% chance of a measurement after each gate) for a 5-qubit system.
"""
function random_ops(N::Int,depth::Int;measure_prob::Float64=0.0,measure_basis::Vector{String}=["MX","MY","MZ"])
opList=Vector{QuantumOps}()#[]#Vector{QuantumOps}(undef,len)

for _=1:depth#layer

    c=1 #runs over qubits
    while c <= N

        if c==N
            s=uppercase(rand(one_qubit_gates))
        else
            s=uppercase(rand(union(one_qubit_gates,two_qubit_gates)))
        end

        p1=round(randn()*pi,digits=2)
        p2=round(randn()*pi,digits=2)
        p3=round(randn()*pi,digits=2)
        
        if s ∈ one_qubit_gates#one qubit
            if s ∈ gates_with_phase
                if s=="U2"
                    push!(opList,Op("$(s)($(p1),$(p2))",c))
                elseif s=="U3"
                    push!(opList,Op("$(s)($(p1),$(p2),$(p3))",c))
                else
                    push!(opList,Op([s,p1],c))
                end
            else
                push!(opList,Op(s,c))
            end

            c=c+1
        elseif s ∈ two_qubit_gates #two qubit
            if s ∈ gates_with_phase
                if s=="FSIM"
                    push!(opList,Op("$(s)($(p1),$(p2))",c,c+1))
                else
                    push!(opList,Op([s,p1],c,c+1))
                end
            else
                push!(opList,Op(s,c,c+1))
            end
            c=c+2
        end

    end


    #measurement layer
    for i=1:N #runs over all sites
        if rand()<measure_prob
            push!(opList,Op(rand(measure_basis),i))
        end
    end

end

return opList
end


"""
    random_ops_len(N::Int,len::Int;measure_prob::Float64=0.0,measure_basis::Vector{String}=["MX","MY","MZ"])

random len not depth
"""
function random_ops_len(N::Int,len::Int;measure_prob::Float64=0.0,measure_basis::Vector{String}=["MX","MY","MZ"])
    opList=Vector{QuantumOps}()#[]#Vector{QuantumOps}(undef,len)
    
    for _=1:len#layer
    
        c=rand(1:N)
    
        if c==N
            s=uppercase(rand(one_qubit_gates))
        else
            s=uppercase(rand(union(one_qubit_gates,two_qubit_gates)))
        end
    
        p1=round(randn()*pi,digits=2)
        p2=round(randn()*pi,digits=2)
        p3=round(randn()*pi,digits=2)
        
        if s ∈ one_qubit_gates#one qubit
            if s ∈ gates_with_phase
                if s=="U2"
                    push!(opList,Op("$(s)($(p1),$(p2))",c))
                elseif s=="U3"
                    push!(opList,Op("$(s)($(p1),$(p2),$(p3))",c))
                else
                    push!(opList,Op([s,p1],c))
                end
            else
                push!(opList,Op(s,c))
            end
    
        elseif s ∈ two_qubit_gates #two qubit
            c2=c==1 ? 1 : rand([1,-1])
            
            if s ∈ gates_with_phase
                if s=="FSIM"
                    push!(opList,Op("$(s)($(p1),$(p2))",c,c+c2))
                else
                    push!(opList,Op([s,p1],c,c+c2))
                end
            else
                push!(opList,Op(s,c,c+c2))
            end
        end
    
        #measurement layer
        for i=1:N #runs over all sites
            if rand()<measure_prob
                push!(opList,Op(rand(measure_basis),i))
            end
        end
    
    end
    
    return opList
end


"""
`random_clifford(N::Int, len::Int; measure_prob::Float64=0.0, measure_basis::Vector{String}=["MX","MY","MZ"]) -> Vector{QuantumOps}`

Generate a random sequence of Clifford gates, with optional mid-circuit measurements.

# Arguments
- `N::Int`: The number of qubits in the system.
- `len::Int`: The length of the sequence of operations to generate.

# Keyword Arguments
- `measure_prob::Float64`: The probability of adding a measurement operation after each gate.
- `measure_basis::Vector{String}`: The basis in which measurements are performed.

# Returns
- `Vector{QuantumOps}`: A vector of randomly chosen quantum operations (`QuantumOps`), each representing a Clifford gate or a measurement operation.

# Description
This function creates a vector of quantum operations, where each operation is either a randomly chosen Clifford gate from the set {"X", "Y", "Z", "H", "S", "CNOT", "SWAP", "CZ"} or a measurement operation, based on `measure_prob`. For two-qubit gates ("CNOT", "SWAP", "CZ"), adjacent qubits (qubit `r` and qubit `r+1`) are selected. For single-qubit gates and measurements, a random qubit `r` is chosen.

# Example
```julia
ops = random_clifford(5, 10; measure_prob=0.2, measure_basis=["MX","MZ"])
```
This example generates a sequence of 10 random Clifford gates and measurements (with a 20% chance of a measurement after each gate) for a 5-qubit system.
"""
function random_clifford(N::Int,len::Int;measure_prob::Float64=0.0,measure_basis::Vector{String}=["MX","MY","MZ"])
    ops=Vector{QuantumOps}()#[]#Vector{QuantumOps}(undef,len)
    for _=1:len
        clifford_rand = rand(["X", "Y", "Z", "H", "S","CNOT", "SWAP", "CZ"])
        c2 = ["CNOT", "SWAP", "CZ"]
        r=rand(1:N-1)

        if clifford_rand ∈ c2
            push!(ops,Op(clifford_rand,r,r+1))
        else
            push!(ops,Op(clifford_rand,r))
        end

        if rand()<measure_prob
            push!(ops,ifOp(rand(measure_basis),rand(1:N)))
        end
    end
    
    return ops
end

"""
`_is_it_measurement(name::String) -> Bool`

Determines if a given gate name represents a measurement operation.

- `name`: Name of the quantum gate.

Returns `true` if the gate is a measurement operation, otherwise `false`.
"""
function _is_it_measurement(name::String)
    uppercase_name=uppercase(name)
    bool=uppercase_name=="MZ" || uppercase_name=="M(Z)" || uppercase_name=="MX" || uppercase_name=="M(X)" || uppercase_name=="MY" || uppercase_name=="M(Y)" || uppercase_name=="MR" || uppercase_name=="M(R)"
    return bool
end

function _get_name_and_arguments(op_name0::String)
    
    op_name=_clean_name(op_name0)

    if BlueTangle._name_with_arg_bool(op_name) 
        op_str = split(op_name0, "(")
        tuple_str = replace(op_str[2], '(' => "", ')' => "")
        elements = string.(split(tuple_str, ","))
        return op_name,"("*op_str[2],eval.(Meta.parse.(elements))
    else
        return op_name
    end

end


_clean_name(op_name::String)=uppercase(string(split(op_name, "(")[1]))

_name_with_arg_bool(name::String)=_clean_name(name) ∈ gates_with_phase
_name_with_two_qubit_gates_bool(name::String)=_clean_name(name) ∈ two_qubit_gates


"""
all gate functions
"""
function gates(op_name::String,param_bool=false)

    _P(λ) = [1 0; 0 exp(im*λ)]

    _RX(θ)=[ #exp(-i θ/2 X)
        cos(θ/2) -im*sin(θ/2);
        -im*sin(θ/2) cos(θ/2)]

    _RY(θ)=[ #exp(-i θ/2 Y)
        cos(θ/2) -sin(θ/2);
        sin(θ/2) cos(θ/2)]

    _RZ(θ)=[ #exp(-i θ/2 Z)
        exp(-im*θ/2) 0;
        0 exp(im*θ/2)]

    _U2(φ, λ)=(1/sqrt(2)) * [1 -exp(im * λ); exp(im * φ) exp(im * (φ + λ))]

    _U3(θ, φ, λ)=[
        cos(θ/2) -exp(im * λ) * sin(θ/2);
        exp(im * φ) * sin(θ/2) exp(im * (φ + λ)) * cos(θ/2)];    

    _CP(lambda)=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 exp(1im*lambda)]

    _GIVENS(theta) = [1 0 0 0; 0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1]

    _FSIM(theta,phi) = [1 0 0 0; 0 cos(theta) -im*sin(theta) 0;0 -im*sin(theta) cos(theta) 0;0 0 0 exp(im*phi)]

    _SWAPA(a)= 0.5 * [2. 0 0 0; 0 1+exp(im*pi*a) 1-exp(im*pi*a) 0; 0 1-exp(im*pi*a) 1+exp(im*pi*a) 0; 0 0 0 2] #SWAP^a

    _RXX(ϕ)=[cos(ϕ) 0 0 -im*sin(ϕ);0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; -im*sin(ϕ) 0 0 cos(ϕ)]

    _RYY(ϕ)=[cos(ϕ) 0 0 im*sin(ϕ);0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; im*sin(ϕ) 0 0 cos(ϕ)]

    _RZZ(ϕ)=[exp(-1im*ϕ/2) 0 0 0;0 exp(1im*ϕ/2) 0 0;0 0 exp(1im*ϕ/2) 0; 0 0 0 exp(-1im*ϕ/2)]

    _RXY(ϕ)=[1 0 0 0;0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; 0 0 0 1]

    local_functions = Dict(
        "P" => _P,
        "RX" => _RX,
        "RY" => _RY,
        "RZ" => _RZ,
        "U2" => _U2,
        "U3" => _U3,
        "CP" => _CP,
        "GIVENS" => _GIVENS,
        "SWAPA" => _SWAPA,
        "FSIM" => _FSIM,
        "RXX" => _RXX,
        "RYY" => _RYY,
        "RZZ" => _RZZ,
        "RXY" => _RXY
    )

    clean_name = _clean_name(op_name)
    uppercase_name = uppercase(op_name)

    if haskey(local_functions, clean_name)

        arg_vec=_get_name_and_arguments(op_name)[3]
        res=local_functions[clean_name](arg_vec...)

        return param_bool==false ?  res : (res,arg_str,arg_vec) 

        #measurement
    elseif uppercase_name=="M(Z)" || uppercase_name=="MZ" || uppercase_name=="M(R)" || uppercase_name=="MR" || uppercase_name=="RES"
        return BlueTangle.gate.I
    elseif uppercase_name=="M(X)" || uppercase_name=="MX"
        return BlueTangle.gate.H
    elseif uppercase_name=="M(Y)" || uppercase_name=="MY"
        return BlueTangle.gate.HSP
    else
        if haskey(BlueTangle.gate, Symbol(uppercase_name))
            return BlueTangle.gate[Symbol(uppercase_name)]
        else
            println("Gate $(op_name) not found")
            return -1
        end
    end

end


# phase gates

export _P,_RX,_RY,_RZ,_U2,_U3,_CP,_GIVENS,_FSIM,_SWAPA,_RXX,_RYY,_RZZ,_RXY

"""
_P(λ) = [1 0; 0 exp(im*λ)]
"""
_P(λ) = [1 0; 0 exp(im*λ)]

_RX(θ)=[ #exp(-i θ/2 X)
    cos(θ/2) -im*sin(θ/2);
    -im*sin(θ/2) cos(θ/2)]

_RY(θ)=[ #exp(-i θ/2 Y)
    cos(θ/2) -sin(θ/2);
    sin(θ/2) cos(θ/2)]

_RZ(θ)=[ #exp(-i θ/2 Z)
    exp(-im*θ/2) 0;
    0 exp(im*θ/2)]

"""
_U2(φ, λ)=(1/sqrt(2)) * [1 -exp(im * λ); exp(im * φ) exp(im * (φ + λ))]
"""
_U2(φ, λ)=(1/sqrt(2)) * [1 -exp(im * λ); exp(im * φ) exp(im * (φ + λ))]

"""
_U3(θ, φ, λ)=[
    cos(θ/2) -exp(im * λ) * sin(θ/2);
    exp(im * φ) * sin(θ/2) exp(im * (φ + λ)) * cos(θ/2)]; 
"""
_U3(θ, φ, λ)=[
    cos(θ/2) -exp(im * λ) * sin(θ/2);
    exp(im * φ) * sin(θ/2) exp(im * (φ + λ)) * cos(θ/2)];    

_CP(lambda)=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 exp(1im*lambda)]

_GIVENS(theta) = [1 0 0 0; 0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1]

_FSIM(theta,phi) = [1 0 0 0; 0 cos(theta) -im*sin(theta) 0;0 -im*sin(theta) cos(theta) 0;0 0 0 exp(im*phi)]

_SWAPA(a)= 0.5 * [2. 0 0 0; 0 1+exp(im*pi*a) 1-exp(im*pi*a) 0; 0 1-exp(im*pi*a) 1+exp(im*pi*a) 0; 0 0 0 2] #SWAP^a

_RXX(ϕ)=[cos(ϕ) 0 0 -im*sin(ϕ);0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; -im*sin(ϕ) 0 0 cos(ϕ)]

_RYY(ϕ)=[cos(ϕ) 0 0 im*sin(ϕ);0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; im*sin(ϕ) 0 0 cos(ϕ)]

_RZZ(ϕ)=[exp(-1im*ϕ/2) 0 0 0;0 exp(1im*ϕ/2) 0 0;0 0 exp(1im*ϕ/2) 0; 0 0 0 exp(-1im*ϕ/2)]

_RXY(ϕ)=[1 0 0 0;0 cos(ϕ) -im*sin(ϕ) 0;0 -im*sin(ϕ) cos(ϕ) 0; 0 0 0 1]