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
- Projectors (`proj0` for |0><0|, `proj1` for |1><1|)
- Controlled gates such as `CX` (CNOT), `CNOT` (an alias for CX), and `CZ`
- The `SWAP`, `iSWAP`, `fSWAP` gate
- The `ECR` gate

Each single-qubit gate is represented as a 2x2 matrix, while multi-qubit gates like `CNOT`, `ECR`, `SYC`, `CZ`, and `SWAP`, `iSWAP`, `fSWAP` are represented as 4x4 matrices.
"""
const gate = (
    I = [1 0; 0 1],
    X  = [0 1; 1 0],
    SX  = (1/2)*[1+1im 1-1im; 1-1im 1+1im],# ==RX(pi/2)*exp(1im*pi/4)
    Y  = [0 -im; im 0],
    Z  = [1 0; 0 -1],
    H  = (1/sqrt(2)) * [1 1; 1 -1],
    S  =[1 0; 0 im],
    T  = [1 0; 0 exp(im * π / 4)],
    HSp = (1/sqrt(2)) * [1+0im  0-1im;1+0im  0+1im],
    proj0=[1 0; 0 0], #|0><0> #projector
    proj1=[0 0; 0 1], #|1><1> #projector
    CX = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], # CNOT
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], # CNOT
    CZ = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1],
    SWAP =[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1],
    iSWAP =[1 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1],
    fSWAP =[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1],
    SYC =[1 0 0 0; 0 0 -im 0; 0 -im 0 0; 0 0 0 exp(-im*pi/6)],#Sycamore
    ECR = (1/sqrt(2)) * [0 + 0im 1 + 0im 0 + 0im 0 + 1im; 1 + 0im 0 + 0im 0 - 1im 0 + 0im; 0 + 0im 0 + 1im 0 + 0im 1 + 0im; 0 - 1im 0 + 0im 1 + 0im 0 + 0im], # (IX-XY)/sqrt(2)
    )

"""
`random_gate_1(N::Int) -> Op`

Generates a random single-qubit gate operation.

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated single-qubit gate.
"""
function random_gate_1(N::Int)## N is total qubit number
    
    lambda=round(rand(),digits=2)
    name=rand(["X","Y","Z","H","S","T","P","RX","RY","RZ"])

    if name in ["P","RX","RY","RZ"]
        name_return="$(name)($(lambda)π)"
    else
        name_return=name
    end

    return Op(name_return,gates1(name,lambda*pi),rand(1:N))
end

"""
`random_gate_2(N::Int) -> Op`

Generates a random two-qubit gate operation.

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated two-qubit gate.
"""
function random_gate_2(N)## N is total qubit number
    
    lambda=round(rand(),digits=2)
    name=rand(["CX","CZ","CP","GIVENS","FSIM","SWAP","iSWAP","fSWAP","SYC","ECR"])

    if name=="CP"
        name_return="$(name)($(lambda)π)"
    elseif name=="GIVENS"
        name_return="$(name)($(lambda)π)"
    elseif name=="FSIM"
        name_return="$(name)($(lambda)π)"
    else
        name_return=name
    end

    r=rand(1:N-1)

    return Op(name_return,gates2(name,lambda*pi),r,r+1)
end

"""
`random_gate(N::Int) -> Op`

Generates a random quantum gate operation (either single or two-qubit).

- `N`: Total number of qubits in the system.

Returns an `Op` object representing the randomly generated quantum gate.
"""
random_gate(N)=rand([random_gate_1(N),random_gate_2(N)])

"""
`random_ops(N::Int, len::Int; measure_prob::Float64=0.0, measure_basis::Vector{String}=["MX", "MY", "MZ"]) -> Vector{QuantumOps}`

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
This function creates a vector of quantum operations, where each operation is either a randomly chosen gate from the set {"X", "Y", "Z", "H", "S", "CX","CZ","CP","GIVENS","FSIM","SWAP","iSWAP","fSWAP","SYC","ECR"} or a measurement operation, based on `measure_prob`. 

# Example
```julia
ops = random_ops(5, 10; measure_prob=0.2, measure_basis=["MX","MZ"])
```
This example generates a sequence of 10 random gates and measurements (with a 20% chance of a measurement after each gate) for a 5-qubit system.
"""
function random_ops(N,len;measure_prob::Float64=0.0,measure_basis::Vector{String}=["MX","MY","MZ"])
    ops=Vector{QuantumOps}()
    for _=1:len
        push!(ops,random_gate(N))
        if rand()<measure_prob
            push!(ops,ifOp(rand(measure_basis),rand(1:N)))
        end
    end
    return ops
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

"""
`_parse_and_call_safe(name::String) -> (String, Array)`

Parses a gate name string with parameters and returns the gate name and parsed parameters.

- `name`: Gate name with potential parameters enclosed in parentheses.

Returns a tuple containing the parsed gate name and an array of parameters.
"""
function _parse_and_call_safe(name::String)
    # Match the function name and the parameter within the parentheses
    match_data = match(r"([A-Za-z]+)\(([^)]+)\)", name)
    
    # Check if the match was successful
    if match_data === nothing || _is_it_measurement(name)
        return name,[]
        print("The gate name format is incorrect: $name")
    end
    
    # Extract the function name and parameter from the regex match
    func_name = string(match_data.captures[1])
    param_str = match_data.captures[2]

    if endswith(uppercase(param_str), "*PI")
        multiplier_str = chop(param_str,tail=3)
        param = parse(ComplexF64, multiplier_str) * pi
    elseif endswith(param_str, "*π")
        multiplier_str = chop(param_str,tail=2)
        param = parse(ComplexF64, multiplier_str) * pi
    elseif (startswith(param_str, "π") && endswith(param_str, "π")) || (startswith(uppercase(param_str), "PI") && endswith(uppercase(param_str), "PI")) #dont change the order
        param = 1.0 * pi
    elseif endswith(param_str, "π")#dont change the order
        multiplier_str = chop(param_str,tail=1)
        param = parse(ComplexF64, multiplier_str) * pi
    elseif endswith(uppercase(param_str), "PI")#dont change the order
        multiplier_str = chop(param_str,tail=2)
        param = parse(ComplexF64, multiplier_str) * pi
    else
        # If there's no "PI" and no "*", parse the number directly
        param = parse(ComplexF64, param_str)
    end

    # Call gates1 function with the parsed function name and parameter
    return func_name,[param]
end


"""
`gates1(op_name::String, param...) -> Matrix`

Constructs a matrix representation of a single-qubit gate based on the gate name and parameters.

- `op_name`: Name of the single-qubit gate.
- `param`: Parameters for the gate, if applicable.

Returns a 2x2 matrix representing the specified single-qubit gate.
"""
function gates1(op_name::String,param...)

    if isempty(param)
        op_name,param=_parse_and_call_safe(op_name)
    end

    uppercase_name = uppercase(op_name)

    P(λ) = [1 0; 0 exp(im*λ)]

    RX(θ)=[ #exp(-i θ/2 X)
        cos(θ/2) -im*sin(θ/2);
        -im*sin(θ/2) cos(θ/2)]

    RY(θ)=[ #exp(-i θ/2 Y)
        cos(θ/2) -sin(θ/2);
        sin(θ/2) cos(θ/2)]

    RZ(θ)=[ #exp(-i θ/2 Z)
        exp(-im*θ/2) 0;
        0 exp(im*θ/2)]
    
    if uppercase_name == "I" || uppercase_name == "ID"
        return gate[:I]
    elseif uppercase_name == "H" || uppercase_name == "HAD"
        return gate[:H]
    elseif uppercase_name == "P" && isempty(param)==false
        return P(param[1])
    elseif uppercase_name == "RX" && isempty(param)==false
        return RX(param[1])
    elseif uppercase_name == "RY" && isempty(param)==false
        return RY(param[1])
    elseif uppercase_name == "RZ" && isempty(param)==false
        return RZ(param[1])

    #measurement
    elseif uppercase_name=="M(Z)" || uppercase_name=="MZ" || uppercase_name=="M(R)" || uppercase_name=="MR"
        return gate.I
    elseif uppercase_name=="M(X)" || uppercase_name=="MX"
        return gate.H
    elseif uppercase_name=="M(Y)" || uppercase_name=="MY"
        return gate.HSp
    
    else
        if haskey(gate, Symbol(uppercase_name))
            return gate[Symbol(uppercase_name)]
        else
            throw("One-qubit gate $(op_name) not found")
        end
    end
end

"""
`gates2(op_name::String, param...) -> Matrix`

Constructs a matrix representation of a two-qubit gate based on the gate name and parameters.

- `op_name`: Name of the two-qubit gate.
- `param`: Parameters for the gate, if applicable.

Returns a 4x4 matrix representing the specified two-qubit gate.
"""
function gates2(op_name::String, param...)

    if isempty(param)
        op_name,param=_parse_and_call_safe(op_name)
    end

    uppercase_name = uppercase(op_name)

    # cx=[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0] # CNOT
    # cz=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]
    cp(lambda)=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 exp(1im*lambda)]
    GIVENS(theta) = [1 0 0 0; 0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1]
    FSIM(theta) = [1 0 0 0; 0 cos(theta) -im*sin(theta) 0;0 -im*sin(theta) cos(theta) 0;0 0 0 1]
    # swap=[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

    if uppercase_name == "CX" || uppercase_name == "CNOT"
        return gate.CX
    elseif uppercase_name == "CZ"
        return gate.CZ # Control-Z gate
    elseif uppercase_name == "SWAP"
        return gate.SWAP
    elseif uppercase_name == "ISWAP"
        return gate.iSWAP
    elseif uppercase_name == "FSWAP"
        return gate.fSWAP
    elseif uppercase_name == "ECR"
        return gate.ECR
    elseif uppercase_name == "SYC"#Sycamore
        return gate.SYC
    elseif uppercase_name == "CP" && length(param)==1
       return cp(param[1])
    elseif uppercase_name == "GIVENS" && length(param)==1
        return GIVENS(param[1])
    elseif uppercase_name == "FSIM" && length(param)==1
        return FSIM(param[1])
    else
        throw("Two-qubit gate $(op_name) not found")
    end
end

# two_qubit_rotation(rθ::Matrix)=kron(rθ, rθ)


"""
`u1(λ)`, `u2(φ, λ)`, `u3(θ, φ, λ)`,

Functions to construct generic single-qubit rotation gates based on Euler angles.

- `u1(λ)`: Rotation about Z-axis.
- `u2(φ, λ)`: Rotation about X+Z axis.
- `u3(θ, φ, λ)`: Generic single-qubit rotation.

Each function returns a 2x2 matrix representing the specified rotation.

Example for u1(λ):
λ=pi -> Z
λ=pi/2 -> S
λ=pi/4 -> T
u1(λ)=U(0, 0, λ)
"""
U1(λ)=[1 0; 0 exp(im * λ)]

"""
`U2(φ, λ)`

Single-qubit rotation about the X+Z axis.

(φ, λ)=(0,pi) -> H
(φ, λ)=(0,0) -> RY(pi/2)
(φ, λ)=(-pi/2,pi/2) -> RX(pi/2)

u2(φ, λ)=U(pi/2, φ, λ)
"""
U2(φ, λ)=(1/sqrt(2)) * [1 -exp(im * λ); exp(im * φ) exp(im * (φ + λ))]

"""

`U3(θ, φ, λ)`

Generic single-qubit rotation gate with 3 Euler angles.

(θ, φ, λ)=(θ, -pi/2, pi/2)=RX(θ)
(θ, φ, λ)=(θ, 0, 0)=RY(θ)
"""
U3(θ, φ, λ)=[
    cos(θ/2) -exp(im * λ) * sin(θ/2);
    exp(im * φ) * sin(θ/2) exp(im * (φ + λ)) * cos(θ/2)
];

# """
# Alias for U3(θ, φ, λ)
# """
# U(θ, φ, λ)=U3(θ, φ, λ)


##========== All gate ==========
