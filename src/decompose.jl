"""
    pauli_decomposition_tensor(n;qubit::Int=0)
Function to generate all tensor products of Pauli matrices for n qubits
"""
function pauli_decomposition_tensor(n;qubit::Int=0,distance::Int=2)

    if distance<=1
        throw(ArgumentError("error with distance."))
    end

    paulis = [sa.sparse(gate.I), sa.sparse(gate.X), sa.sparse(gate.Y), sa.sparse(gate.Z)]
    indices = Iterators.product(ntuple(_ -> 1:4, n)...)

    if qubit==0

        tensor_products = []
        for index in indices
            pauli_matrices = [paulis[i] for i in index]
            push!(tensor_products, foldl(kron,(pauli_matrices)))
        end

    else
        
        products = []
        for index in indices
            pauli_matrices = [paulis[i] for i in index]
            push!(products, pauli_matrices)
        end

        tensor_products = []
        for product=products

            if qubit>1
                for _=1:qubit-1
                    insert!(product,1,paulis[1])
                end
            end

            for _=1:distance-1
                insert!(product,qubit+1,paulis[1])
            end

            push!(tensor_products, foldl(kron,product))
        end

    end

    return tensor_products

end

"""
    pauli_decomposition_names(n::Int)
Function to generate names of all tensor products of Pauli matrices for n qubits
"""
function pauli_decomposition_names(n::Int)
    pauli_names=["I", "x", "y", "z"]
    indices = Iterators.product(ntuple(_ -> 1:4, n)...)

    tensor_products = []
    for index in indices
            push!(tensor_products, join([pauli_names[i] for i in index], "*"))
    end

    return tensor_products
end

"""
    pauli_decomposition(A::AbstractMatrixS)
"""
function pauli_decomposition(A::AbstractMatrixS)
    n = Int(log2(size(A, 1)))
    pauli_tensors = pauli_decomposition_tensor(n)
    coefficients = [la.tr(A * P) / 2^n for P in pauli_tensors]
    return sa.sparse(coefficients), pauli_tensors
end

"""
pauli_reconstruction(coefficients::AbstractVectorS, pauli_tensors::Vector)
# Function to reconstruct the matrix from Pauli decomposition

pauli_reconstruction(coefficients::AbstractVectorS,qubit::Int)
    # nonlocal reconstruction #add id on qubit+1
"""
function pauli_reconstruction(coefficients::AbstractVectorS, pauli_tensors::Vector)
    sum(coefficients[i] * pauli_tensors[i] for i in 1:length(coefficients))
end

#nonlocal reconstruction
function pauli_reconstruction(coefficients::AbstractVectorS,qubit::Int;distance::Int=2)
    n=Int(log(4,length(coefficients)))
    pauli_nonlocal_tensors=pauli_decomposition_tensor(n;qubit=qubit,distance=distance)
    sum(coefficients[i] * pauli_nonlocal_tensors[i] for i in 1:length(coefficients))
end


##
##========== decomposition ==========

"""
zyz_decomposition(U::AbstractMatrixS)

U = gate.H
α, β, γ, δ = zyz_decomposition(U)
U2 = exp(im*α) * _RZ(β) * _RY(γ) * _RZ(δ)
isapprox(la.norm(U - U2),0,atol=1e-10)
"""
function zyz_decomposition(U::AbstractMatrixS)

    # https://quantumcomputing.stackexchange.com/a/16263

    if !isunitary(U)
        throw("U is not unitary!")
    end

    a, b = U[1, 1], U[1, 2]
    c, d = U[2, 1], U[2, 2]

    # Determine γ
    if a != 0
        γ = 2 * atan(abs(b) / abs(a))
    else
        γ = π
    end

    # Determine β and δ
    if γ == 0
        β_plus_δ = angle(d) - angle(a)
        β = 0
        δ = β_plus_δ
    elseif γ == π
        β_minus_δ = angle(-b) - angle(c)
        β = 0
        δ = -β_minus_δ
    else
        β = angle(c) - angle(a)
        δ = angle(-b) - angle(a)
    end

    # Determine α
    if a != 0
        α = angle(a) + (β + δ) / 2
    else
        α = angle(c) + (-β + δ) / 2
    end

    return (α, β, γ, δ) # = exp(im*α)*Rz(β) * Ry(γ) * Rz(δ)
end


"""
kronecker_decomposition(C::AbstractMatrixS)

C = kron(gates("RZ(.1pi)"), gates("RY(.3)"))
A, B = nearest_kronecker_product(C)
isapprox(la.norm(kron(A, B) - C),0,atol=1e-10)
"""
function kronecker_decomposition(C::AbstractMatrixS)

    #after Eq. 87 from https://github.com/gecrooks/on_gates

    # Reshape C from a 4x4 matrix into a 2x2x2x2 tensor
    C = reshape(C, 2, 2, 2, 2)
    
    # Permute the dimensions of C to match the outer product structure
    C = permutedims(C, (1, 3, 2, 4))
    
    # Flatten C back into a 4x4 matrix
    C = reshape(C, 4, 4)

    U, S, V = la.svd(C)

    A = sqrt(S[1]) * reshape(V[:, 1]', 2, 2)  # Transpose V with ' for correct dimensions in Julia
    B = sqrt(S[1]) * reshape(U[:, 1], 2, 2)

    return A,B
end

##========== decomposition ==========