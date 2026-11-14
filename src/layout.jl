# using SparseArrays

function illustrate_physical_qubits(qubit_grid::sa.SparseMatrixCSC{Int, Int}, qubit_num_to_pos)
    num_rows, num_cols = size(qubit_grid)
    illustrated_grid = zeros(Int, num_rows, num_cols) # Initialize a grid for illustration

    for (qubit_num, pos) in qubit_num_to_pos
        illustrated_grid[pos] = qubit_num # Place qubit numbers in their positions
    end

    # println("Physical Qubits in the Grid:")
    return sa.sparse(illustrated_grid)
end

# Function to enumerate physical qubits and create mappings
function enumerate_qubits(qubit_grid::sa.SparseMatrixCSC{Int, Int})
    qubit_number = 0
    qubit_num_to_pos = Dict{Int, CartesianIndex{2}}()
    pos_to_qubit_num = Dict{CartesianIndex{2}, Int}()

    for pos in findall(qubit_grid .> 0)
        qubit_number += 1
        qubit_num_to_pos[qubit_number] = pos
        pos_to_qubit_num[pos] = qubit_number
    end

    return qubit_num_to_pos, pos_to_qubit_num
end

# Function to generate a neighbors dictionary for each physical qubit
function generate_neighbors_dict(qubit_num_to_pos, pos_to_qubit_num)
    neighbors_dict = Dict{Int, Set{Int}}()
    for (qubit_num, pos) in qubit_num_to_pos
        neighbors = Set{Int}()
        
        # Check potential neighbors by applying Cartesian offsets
        for delta in [CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(0, 1)]
            neighbor_pos = pos + delta
            if neighbor_pos in keys(pos_to_qubit_num)
                neighbor_qubit_num = pos_to_qubit_num[neighbor_pos]
                push!(neighbors, neighbor_qubit_num)
            end
        end
        
        neighbors_dict[qubit_num] = neighbors
    end

    return neighbors_dict
end

# Function to check if two qubits are neighbors
_are_neighbors(qubit_num1::Int, qubit_num2::Int, neighbors_dict::Dict)=qubit_num1 ∈ neighbors_dict[qubit_num2]

function _average_connectivity(qubit_grid::sa.SparseMatrixCSC{Int, Int})
    # Enumerate qubits and generate mappings
    qubit_num_to_pos, pos_to_qubit_num = enumerate_qubits(qubit_grid)
    
    # Generate neighbors dictionary
    neighbors_dict = generate_neighbors_dict(qubit_num_to_pos, pos_to_qubit_num)
    
    # Calculate the total number of connections and the average connectivity
    total_connections = sum(length(neighbors) for neighbors in values(neighbors_dict))
    num_qubits = length(qubit_num_to_pos)
    average_connectivity = total_connections / num_qubits
    
    return average_connectivity
end

function _make_swaps(op::QuantumOps, neighbors::Dict; noisy_swap::Bool=false, fswap::Bool=false)
    op isa OpF && return [op]
    if op isa ifOp
        branches = [Op[routed for branch_op in branch for routed in _make_swaps(branch_op, neighbors; noisy_swap=noisy_swap, fswap=fswap)] for branch in op.if01]
        return [ifOp(op.name, op.qubit, branches[1], branches[2])]
    end
    if op isa Op && op.q == 1 && op.control > 0
        op.qubit in neighbors[op.control] && return [op]
        path = _find_shortest_path(neighbors, op.control, op.qubit)
        isempty(path) && throw(ArgumentError("the control and target are disconnected in this layout"))
        swap_name = fswap ? "FSWAP" : "SWAP"
        swaps = [Op(swap_name, path[i], path[i+1]; noisy=noisy_swap) for i in 1:length(path)-2]
        routed = Op(op.name, op.mat, op.qubit; noisy=op.noisy, control=path[end-1])
        return [swaps; routed; reverse(swaps)]
    end
    (op.q == 1 || op.q == 3 || op.qubit in neighbors[op.target_qubit]) && return [op]
    path = _find_shortest_path(neighbors, op.qubit, op.target_qubit)
    isempty(path) && throw(ArgumentError("the operated qubits are disconnected in this layout"))
    swap_name = fswap ? "FSWAP" : "SWAP"
    swaps = [Op(swap_name, path[i], path[i+1]; noisy=noisy_swap) for i in 1:length(path)-2]
    control = op.control
    for swap in swaps
        control = control == swap.qubit ? swap.target_qubit : control == swap.target_qubit ? swap.qubit : control
    end
    # Preserve the actual matrix/function/channel, including custom gates and
    # operand order; routing is a permutation followed by its inverse.
    routed = if op isa OpQC
        OpQC(op.name, op.kraus, path[end-1], op.target_qubit; type=op.type)
    else
        Op(op.name, op.mat, path[end-1], op.target_qubit; noisy=op.noisy, control=control)
    end
    return [swaps; routed; reverse(swaps)]
end


function _find_shortest_path(neighbors_dict, start_qubit, end_qubit)
    visited = Set([start_qubit])
    queue = [(start_qubit, [start_qubit])]  # (current_qubit, path_so_far)
    
    next = 1
    while next <= length(queue)
        current_qubit, path = queue[next]
        next += 1
        if current_qubit == end_qubit
            return path
        end
        for neighbor in neighbors_dict[current_qubit]
            if !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, (neighbor, [path; neighbor]))
            end
        end
    end
    
    return []  # Return an empty path if no path is found
end
