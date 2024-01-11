"""
`noise_model1(model::String, p::Float64)` and `noise_model2(model::String, p::Float64)`

Wrappers for creating single-qubit and two-qubit noise models, respectively.

- `model`: Name of the noise model.
- `p`: Parameter for the noise model.

Returns a noise model for either one or two qubits, depending on the function used.
"""
noise_model1(model::String, p::Float64)=noise_model(model, p; two_qubit=false)
noise_model2(model::String, p::Float64)=noise_model(model, p; two_qubit=true)

"""
`apply_noise(ops::Vector{T},noise12::Tuple) where T <: QuantumOps -> Vector{QuantumOps}`

Applies specified noise models to a vector of quantum operations.

- `ops`: Vector of QuantumOps representing the operations.
- `noise12`: A tuple containing the single-qubit and two-qubit noise models.

Returns a new vector of QuantumOps with noise models applied.
"""
function apply_noise(ops::Vector{T},noise12::Tuple{QuantumChannel, QuantumChannel}) where T <: QuantumOps

    noise1,noise2=noise12

    if noise1.q==1 && noise2.q==2
        new_ops=Vector{QuantumOps}()
        for op in ops
            
            if op.noise==false #do not override noise!
                if op.q==1
                    push!(new_ops,Op(op.name,op.gate,op.qubit,noise1))
                elseif op.q==2
                    push!(new_ops,Op(op.name,op.gate,op.qubit,op.target_qubit,noise2))
                else
                    throw("only works for 1 or 2 qubit operators")
                end
            end

        end

    else
        throw("single qubit or two-qubit noise is not correct")
    end

    return new_ops
end

"""
    noise_model(model::String, p::Union{Float64,ComplexF64}; two_qubit::Bool=false) -> Array

Generates a quantum noise model based on the specified type and parameters.

# Arguments
- `model`: Name of the noise model. Available models:
  - `amplitude_damping`: Simulates the loss of energy from a quantum system. The parameter `γ` (given by `p`) represents the probability of the qubit losing its excitation.
  - `phase_damping`: Represents the loss of quantum information without loss of energy. The parameter `γ` quantifies the likelihood of a phase shift occurring.
  - `phase_flip`: Introduces a phase error (Z error). The qubit's phase is flipped with probability `p`.
  - `bit_flip`: Causes a flip in the qubit's state (X error) with probability `p`.
  - `bit_phase_flip`: Combines bit flip and phase flip errors, effectively applying a Y error with probability `p`.
  - `depolarizing`: A general error model where any of the Pauli operations (X, Y, Z) can be applied with equal probability. The parameter `p` represents the overall error probability.
  - `rot_X`, `rot_Y`, `rot_Z`, `rot_P`: Represent coherent errors (incorrect rotations) around the X, Y, Z, and phase axes, respectively. The parameter `p` specifies the rotation error magnitude.
- `p`: Parameter for the noise model, representing probabilities or magnitudes of errors.
- `two_qubit`: Boolean flag to indicate if the noise model is for two qubits. When set to `true`, the model is extended to two-qubit operations.

# Returns
An array of operations representing the noise model.

# Examples
To create a single-qubit phase damping noise model with a damping parameter of 0.05:
```julia
noise = Noise1("phase_damping", 0.01)
```

For a two-qubit depolarizing error model with an error probability of 0.1:

```julia
noise = Noise2("phase_damping", 0.1)
```
"""
function noise_model(model::String, p::Union{Float64,ComplexF64}; two_qubit=false)

    id = gate.I
    X = gate.X
    Y = gate.Y
    Z = gate.Z

    rz(θ)=[exp(-im*θ/2) 0; 0 exp(im*θ/2)];
    ry(θ)=[cos(θ/2) -sin(θ/2);sin(θ/2) cos(θ/2)]
    rx(θ)=[cos(θ/2) -im*sin(θ/2);-im*sin(θ/2) cos(θ/2)]
    Ph(λ) = [1 0; 0 exp(im*λ)]

    model=lowercase(model)

    if model == "amplitude_damping"
        γ = p #note that gamma is not same as probability
        E₀ = [1.0 0; 0 sqrt(1 - γ)]
        E₁ = [0 sqrt(γ); 0 0]
        ops=[E₀, E₁]
        
    elseif model == "phase_damping"
        γ = p #note that gamma is not same as probability
        E₀ = [1.0 0; 0 sqrt(1 - γ)]
        E₁ = [0 0; 0 sqrt(γ)]
        ops=[E₀, E₁]
        
    elseif model == "phase_flip" || model == "bit_flip" || model == "bit_phase_flip"
        E₀ = sqrt(1 - p) * id
        E₁ = sqrt(p) * (model == :phase_flip ? Z : (model == :bit_flip ? X : Y))
        ops=[E₀, E₁]
        
    elseif model == "depolarizing_amp"#|| model == "depolarizing_qiskit"
        E₀ = sqrt(1 - 3p/4) * id
        E₁ = sqrt(p/4) * X
        E₂ = sqrt(p/4) * Y
        E₃ = sqrt(p/4) * Z
        ops=[E₀, E₁, E₂, E₃]

    elseif model == "depolarizing"
        E₀ = sqrt(1 - p) * id
        E₁ = sqrt(p/3) * X
        E₂ = sqrt(p/3) * Y
        E₃ = sqrt(p/3) * Z
        ops=[E₀, E₁, E₂, E₃]

    elseif model == "rot_z" #coherent error (incorrect rotation)
        ops=[gate.I,rz(pi*p)]/sqrt(2)
        
    elseif model == "rot_y" #coherent error (incorrect rotation)
        ops=[gate.I,ry(pi*p)]/sqrt(2)

    elseif model == "rot_x" #coherent error (incorrect rotation)
        ops=[gate.I,rx(pi*p)]/sqrt(2)
        
    elseif model == "rot_p" #coherent error (incorrect rotation)
        ops=[gate.I,Ph(pi*p)]/sqrt(2)

    elseif model == "rot_err"
        ops=[gate.I,rx(pi*p),ry(pi*p),rz(pi*p)]/sqrt(4)

    #measurements as a quantum channel for classical_shadow
    elseif _is_it_measurement(model)
        ops=[gate.proj0, gate.proj1]#these are projectors
    else
        throw(ArgumentError("Unknown quantum error model"))
    end

    if two_qubit==false
        return ops
    else
        # if model=="depolarizing_qiskit"
                # ops2=qiskit_depolarizing_two_qubit(p)
        # else
                ops2=[kron(Ki, Kj) for Ki in ops for Kj in ops]
        # end

        return ops2
    end
end

##
##========== Twirl ==========

"""
`apply_twirl(op::Op) -> Vector{QuantumOps}`

Applies twirling to a single quantum operation. This is used to randomize errors in quantum operations.

- `op`: A QuantumOps object representing the quantum operation.

Returns a vector of QuantumOps representing the twirled operation.
"""
apply_twirl(op::Op)=apply_twirl(op,false)

"""
`apply_twirl(op::Op, twirling_noise::Union{QuantumChannel,Bool}) -> Vector{QuantumOps}`

Applies twirling with specified noise to a single quantum operation.

- `op`: A QuantumOps object representing the quantum operation.
- `twirling_noise`: Noise model to be applied during twirling.

Returns a vector of QuantumOps representing the twirled operation with noise.
"""
function apply_twirl(op::Op,twirling_noise::Union{QuantumChannel,Bool})

    if op.q==1
        return [op]
        # throw("error! works only for two qubit gates")
    end

    if op.noise==false
        println("no error in op why twirling? -- $(op.name)[$(op.qubit),$(op.target_qubit)]!")
        return [op]
    end

    qubit=op.qubit
    target_qubit=op.target_qubit
    least=min(qubit,target_qubit)
    gate=op.gate

    id = [1 0; 0 1]  # Identity
    X = [0 1; 1 0]  # Pauli-X gate
    Y = [0 -im; im 0]  # Pauli-Y gate
    Z = [1 0; 0 -1]  # Pauli-Z gate

    # Map from symbols to matrices
    pauli_dict = Dict("I" => id, "X" => X, "Y" => Y, "Z" => Z)

    list_of_two_qubits = []
    for op1 in keys(pauli_dict), op2 in keys(pauli_dict)
        op1_matrix = pauli_dict[op1]
        op2_matrix = pauli_dict[op2]
        push!(list_of_two_qubits, (kron(op1_matrix, op2_matrix), (op1,op2)))
    end

    twirl_set=Set()
    for (p1_k,p1_s) in list_of_two_qubits
        for (p2_k,p2_s) in list_of_two_qubits

            if p1_k*gate*p2_k==gate ||  p1_k*gate*p2_k==-1*gate ||  p1_k*gate*p2_k==1im*gate ||  p1_k*gate*p2_k==-1im*gate

                before=[Op(p1_s[1],least,twirling_noise),Op(p1_s[2],least+1,twirling_noise)]
                after=[Op(p2_s[1],least,twirling_noise),Op(p2_s[2],least+1,twirling_noise)]

                push!(twirl_set,(before,op,after))
            end

        end
    end
    return vcat(rand(twirl_set)...)

end

"""
`apply_twirl(ops::Vector{QuantumOps}, twirling_noise::Union{QuantumChannel,Bool}) -> Vector{QuantumOps}`

Applies twirling with specified noise to a vector of quantum operations.

- `ops`: Vector of QuantumOps representing the operations.
- `twirling_noise`: Noise model to be applied during twirling.

Returns a vector of QuantumOps with each operation twirled with the specified noise.
"""
function apply_twirl(ops::Vector{QuantumOps},twirling_noise::Union{QuantumChannel,Bool})
    list_of_all_twirl=Vector{QuantumOps}()
    for op in ops
        append!(list_of_all_twirl,apply_twirl(op,twirling_noise))
    end
    return list_of_all_twirl
end


##========== Twirl ==========


# make correlated noise model