"""
`noise_model1(model::String, p::Float64)` and `noise_model2(model::String, p::Float64)`

Wrappers for creating single-qubit and two-qubit noise models, respectively.

- `model`: Name of the noise model.
- `p`: Parameter for the noise model.

Returns a noise model for either one or two qubits, depending on the function used.
"""
noise_model1(model::String, p::Float64)=noise_model(model, p; two_qubit=false)


"""
noise_model2(model::String, p::Float64)=noise_model(model, p; two_qubit=true)
"""
noise_model2(model::String, p::Float64)=noise_model(model, p; two_qubit=true)

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
  - `rot_X`, `rot_Y`, `rot_Z`, `rot_P`: Represent coherent errors (incorrect rotations) around the X, Y, Z, and phase axes, respectively. The parameter `p` specifies the rotation error magnitude. Note that `p` is multiplied by π.
  - `rot_xyz`: Rotation error in X,Y,Z with equal probability.
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
        E₁ = sqrt(p) * (model == :phase_flip ? Z : (model == "bit_flip" ? X : Y))
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
        ops=[rz(p)]
        
    elseif model == "rot_y" #coherent error (incorrect rotation)
        ops=[ry(p)]

    elseif model == "rot_x" #coherent error (incorrect rotation)
        ops=[rx(p)]
        
    elseif model == "rot_p" #coherent error (incorrect rotation)
        ops=[Ph(p)]

    elseif model == "rot_xyz"
        ops=[rx(p),ry(p),rz(p)]/sqrt(3)

    #measurements as a quantum channel for classical_shadow
    elseif _is_it_measurement(model)
        ops=[gate.P0, gate.P1]#these are projectors
        
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


##========== Twirl ==========

function _twirl_set(name::String,qubit::Int,target_qubit::Int)

    list_ecr=[
    ["Y", "Y", "Y", "Y"],
    ["Y", "Z", "Z", "I"],
    ["Y", "I", "Z", "Z"],
    ["Y", "X", "Y", "X"],
    ["Z", "Y", "Z", "Y"],
    ["Z", "Z", "Y", "I"],
    ["Z", "I", "Y", "Z"],
    ["Z", "X", "Z", "X"],
    ["I", "Y", "X", "X"],
    ["I", "Z", "I", "Z"],
    ["I", "I", "I", "I"],
    ["I", "X", "X", "Y"],
    ["X", "Y", "I", "X"],
    ["X", "Z", "X", "Z"],
    ["X", "I", "X", "I"],
    ["X", "X", "I", "Y"]
    ]

    list_ecr_rev=[
    ["Y", "Y", "Y", "Y"],
    ["Y", "Z", "Y", "Z"],
    ["Y", "I", "X", "X"],
    ["Y", "X", "X", "I"],
    ["Z", "Y", "I", "Z"],
    ["Z", "Z", "I", "Y"],
    ["Z", "I", "Z", "I"],
    ["Z", "X", "Z", "X"],
    ["I", "Y", "Z", "Z"],
    ["I", "Z", "Z", "Y"],
    ["I", "I", "I", "I"],
    ["I", "X", "I", "X"],
    ["X", "Y", "X", "Y"],
    ["X", "Z", "X", "Z"],
    ["X", "I", "Y", "X"],
    ["X", "X", "Y", "I"]
    ]

    list_cnot=[
     ["Y", "Y", "X", "Z"],
     ["Y", "Z", "X", "Y"],
     ["Y", "I", "Y", "X"],
     ["Y", "X", "Y", "I"],
     ["Z", "Y", "I", "Y"],
     ["Z", "Z", "I", "Z"],
     ["Z", "I", "Z", "I"],
     ["Z", "X", "Z", "X"],
     ["I", "Y", "Z", "Y"],
     ["I", "Z", "Z", "Z"],
     ["I", "I", "I", "I"],
     ["I", "X", "I", "X"],
     ["X", "Y", "Y", "Z"],
     ["X", "Z", "Y", "Y"],
     ["X", "I", "X", "X"],
     ["X", "X", "X", "I"]]
    
    list_cnot_rev=[
     ["Y", "Y", "Z", "X"],
     ["Y", "Z", "Y", "I"],
     ["Y", "I", "Y", "Z"],
     ["Y", "X", "Z", "Y"],
     ["Z", "Y", "Y", "X"],
     ["Z", "Z", "Z", "I"],
     ["Z", "I", "Z", "Z"],
     ["Z", "X", "Y", "Y"],
     ["I", "Y", "X", "Y"],
     ["I", "Z", "I", "Z"],
     ["I", "I", "I", "I"],
     ["I", "X", "X", "X"],
     ["X", "Y", "I", "Y"],
     ["X", "Z", "X", "Z"],
     ["X", "I", "X", "I"],
     ["X", "X", "I", "X"]
    ]
    
    if qubit<target_qubit
        return uppercase(name)=="ECR" ? rand(list_ecr) : rand(list_cnot)
    else
        return uppercase(name)=="ECR" ? rand(list_ecr_rev) : rand(list_cnot_rev)
    end
    
end

"""
`apply_twirl(op::QuantumOps) -> Vector{QuantumOps}`

Applies twirling with specified noise to a single quantum operation.

- `op`: A QuantumOps object representing the quantum operation.

Returns a vector of QuantumOps representing the twirled operation with noise.
"""
function apply_twirl(op::QuantumOps)

    if op.q==2 && (uppercase(op.name) == "CNOT" || uppercase(op.name) =="CX" || uppercase(op.name) == "ECR")

        local_ops=Vector{QuantumOps}()
        
        least=min(op.qubit,op.target_qubit)

        twirl_list=_twirl_set(op.name,op.qubit,op.target_qubit)

        for (i,set_op)=enumerate(twirl_list)
            final_qubit=isodd(i) ? least : least+1
            if i==3
                push!(local_ops,op)
            end
            push!(local_ops,Op(set_op,final_qubit))
        end
        
        return local_ops

    else

        return [op]
    end

end

    