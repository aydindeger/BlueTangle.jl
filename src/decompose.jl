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
    pauli_decomposition(A::Union{AbstractMatrix,sa.SparseMatrixCSC})
"""
function pauli_decomposition(A::Union{AbstractMatrix,sa.SparseMatrixCSC})
    n = Int(log2(size(A, 1)))
    pauli_tensors = pauli_decomposition_tensor(n)
    coefficients = [la.tr(A * P) / 2^n for P in pauli_tensors]
    return sa.sparse(coefficients), pauli_tensors
end

"""
pauli_reconstruction(coefficients::sa.SparseVector, pauli_tensors::Vector)
# Function to reconstruct the matrix from Pauli decomposition

pauli_reconstruction(coefficients::sa.SparseVector,qubit::Int)
    # nonlocal reconstruction #add id on qubit+1
"""
function pauli_reconstruction(coefficients::sa.SparseVector, pauli_tensors::Vector)
    sum(coefficients[i] * pauli_tensors[i] for i in 1:length(coefficients))
end

#nonlocal reconstruction
function pauli_reconstruction(coefficients::sa.SparseVector,qubit::Int;distance::Int=2)
    n=Int(log(4,length(coefficients)))
    pauli_nonlocal_tensors=pauli_decomposition_tensor(n;qubit=qubit,distance=distance)
    sum(coefficients[i] * pauli_nonlocal_tensors[i] for i in 1:length(coefficients))
end

