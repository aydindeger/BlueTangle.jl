"""
    trotter_ising(t::Float64, N::Int, J::Float64, h::Float64; dt=0.1) -> Vector{QuantumOps}

Simulates the Ising model Hamiltonian using Trotterization.

- `total_time`: The total time for the simulation.
- `N`: Number of qubits in the Ising chain.
- `J`: Coupling constant for the interaction between adjacent qubits.
- `h`: Magnetic field strength applied in the x-direction.
- `dt` (optional): Time step for the Trotterization. Default is 0.1.

Performs Trotterization of the Ising model Hamiltonian over `N` qubits for a total time `total_time`, given the coupling constant `J` and magnetic field strength `h`. The function constructs a series of quantum gates that approximate the evolution under the Ising Hamiltonian.

Returns a vector of `QuantumOps`, representing the sequence of operations for the Trotterized Ising model simulation.
"""
function trotter_ising(N::Int,total_time::Float64,J::Float64,h::Float64;dt=0.1)
    trotter_step=round(Int,total_time/dt)
    x_angle=round(-2*dt*h,sigdigits=6)
    z_angle=round(-2*dt*J,sigdigits=6)
    final_list=vcat([Op("RX($(x_angle))",i) for i=1:N],[[Op("CNOT",i,i+1),Op("RZ($(z_angle))",i+1),Op("CNOT",i,i+1)] for i=1:N-1]...)
    
    return vcat(fill(final_list,trotter_step)...)#list of operators
end

"""
    hamiltonian_exp(N::Int, total_time::Float64, string_of_ops::Vector; dt=0.1) -> Vector{QuantumOps}

Expands a given Hamiltonian expressed as a string of operations into a sequence of quantum gates using Trotterization.

- `N`: Number of qubits involved in the simulation.
- `total_time`: The total simulation time over which the Hamiltonian is to be applied.
- `string_of_ops`: A vector where odd indices contain the coupling strengths and even indices contain comma-separated strings representing the operators (X, Y, Z) applied to consecutive qubits.
- `dt` (optional): The time step for Trotterization, with a default value of 0.1.

The function `hamiltonian_exp` parses the `string_of_ops` to construct a sequence of operations based on the specified operators and their coupling strengths. For each term in the Hamiltonian:

This method returns a vector of `QuantumOps`, each representing a quantum operation to be applied sequentially to simulate the Hamiltonian over the specified time `total_time`.

### Example
# Define a Hamiltonian for a 3-qubit system with mixed interactions over 2 seconds
N = 3
total_time = 2.0
string_of_ops = [1.0, "X,Y", 0.5, "Y,Z"]
ops = hamiltonian_exp(N, total_time, string_of_ops)

# ops will contain a sequence of quantum gates to apply.
"""
function hamiltonian_exp(N::Int, total_time::Float64, string_of_ops::Vector; dt=0.01, full=true)
    couplings = string_of_ops[1:2:end]
    ops = string_of_ops[2:2:end]
    
    term_ops = [split(op, ",") for op in ops]
    len_op = [length(split(op, ",")) for op in ops]
    
    trotter_step = round(Int, total_time/dt)
    
    all_ops = Vector{QuantumOps}()
    
    for term in 1:length(ops)
        terms = term_ops[term]
        j = len_op[term] - 1
        angle = round(-2 * dt * couplings[term], sigdigits=6)
        
        for i in 1:N-j
            _apply_term(all_ops, terms, angle, collect(i:i+j))
        end
    end
    
    return full ? vcat(fill(all_ops, trotter_step)...) : all_ops
end

function _apply_term(all_ops::Vector{QuantumOps}, terms::Vector{SubString{String}}, angle::Float64, qubits::Vector{Int})

    for (en, qubit) in enumerate(qubits)
        if terms[en] == "X"
            push!(all_ops, Op("H", qubit))
        elseif terms[en] == "Y"
            push!(all_ops, Op("HY", qubit))
        end
    end
    
    # Apply CNOT ladder
    for i in 1:length(qubits)-1
        push!(all_ops, Op("CNOT", qubits[i], qubits[i+1]))
    end
    
    # Apply RZ rotation on the last qubit
    push!(all_ops, Op("RZ($angle)", qubits[end]))
    
    # Apply reverse CNOT ladder
    for i in length(qubits)-1:-1:1
        push!(all_ops, Op("CNOT", qubits[i], qubits[i+1]))
    end
    
    # Apply final basis change
    for (en, qubit) in enumerate(qubits)
        if terms[en] == "X"
            push!(all_ops, Op("H", qubit))
        elseif terms[en] == "Y"
            push!(all_ops, Op("HY", qubit)')
        end
    end
end

function _get_qubit_index(rows::Int, cols::Int, row::Int, col::Int, style::Symbol)
    if style == :rowwise
        return (row - 1) * cols + col
    elseif style == :colwise
        return (col - 1) * rows + row
    else
        throw(ArgumentError("Unsupported enumeration style"))
    end
end

function _create_lattice_map(rows::Int, cols::Int, enumeration_style::Symbol=:rowwise)
    if enumeration_style == :rowwise
        return [[(r-1)*cols + c for c in 1:cols] for r in 1:rows]
    elseif enumeration_style == :colwise
        return [[r + (c-1)*rows for r in 1:rows] for c in 1:cols]
    else
        throw(ArgumentError("Unsupported enumeration style. Use :rowwise or :colwise."))
    end
end

"""
    hamiltonian_exp(rows_cols::Union{Tuple{Int64, Int64}, Vector{Int64}}, total_time::Float64, string_of_ops::Vector;
                    dt=0.01, full=true, enumeration_style::Symbol=:rowwise) -> Vector{QuantumOps}

Generate a sequence of quantum operations to simulate a 2D lattice Hamiltonian using Trotterization.

# Arguments
- `rows_cols`: A tuple or vector specifying the dimensions of the 2D lattice (rows, columns).
- `total_time`: The total simulation time.
- `string_of_ops`: A vector alternating between coupling strengths and operator strings. 
  Each operator string is a comma-separated list of Pauli operators (X, Y, Z) representing 
  the terms in the Hamiltonian.

# Keyword Arguments
- `dt`: Time step for Trotterization (default: 0.01).
- `full`: If true, repeats the operation sequence for each Trotter step (default: true).
- `enumeration_style`: Specifies how qubits are indexed in the 2D lattice. 
  Options are `:rowwise` (default) or `:colwise`.

# Returns
A vector of `QuantumOps` representing the sequence of quantum operations for the simulation.

# Details
This function implements a Trotterized evolution of a 2D lattice Hamiltonian. It supports
arbitrary Pauli string operators and handles both horizontal and vertical nearest-neighbor 
interactions in the lattice.

The Hamiltonian terms are applied to all relevant qubit pairs in the lattice. For operators
involving two qubits (e.g., "X,Y"), both horizontal and vertical applications are considered.
Single-qubit operators are applied to each qubit individually.

The function uses the `_apply_term` helper function to generate the specific quantum operations
for each term in the Hamiltonian.

# Example
```julia
rows, cols = 3, 4
total_time = 2.0
string_of_ops = [1.0, "X,Y", 0.5, "Z"]
ops = hamiltonian_exp((rows, cols), total_time, string_of_ops)
```
"""
function hamiltonian_exp(rows_cols::Union{Tuple{Int64, Int64},Vector{Int64}}, total_time::Float64, string_of_ops::Vector; dt=0.01, full=true, enumeration_style::Symbol=:rowwise)
    rows, cols = rows_cols
    N = rows * cols
    
    couplings = string_of_ops[1:2:end]
    ops = string_of_ops[2:2:end]
    
    term_ops = [split(op, ",") for op in ops]
    len_op = [length(split(op, ",")) for op in ops]
    
    trotter_step = round(Int, total_time/dt)
    
    all_ops = Vector{QuantumOps}()
    
    for term in 1:length(ops)
        terms = term_ops[term]
        j = len_op[term] - 1
        angle = round(-2 * dt * couplings[term], sigdigits=6)
        
        # horizontal connections
        for row in 1:rows
            for col in 1:cols-j
                qubits = [_get_qubit_index(rows, cols, row, c, enumeration_style) for c in col:col+j]
                _apply_term(all_ops, terms, angle, qubits)
            end
        end
        
        # vertical connections if j == 1 (nearest neighbor)
        if j == 1
            for row in 1:rows-1
                for col in 1:cols
                    qubits = [_get_qubit_index(rows, cols, r, col, enumeration_style) for r in row:row+1]
                    _apply_term(all_ops, terms, angle, qubits)
                end
            end
        end
    end
    
    return full ? vcat(fill(all_ops, trotter_step)...) : all_ops
end

"""
test this then delete
"""
function hamiltonian_exp_test(rows_cols::Union{Tuple{Int64, Int64},Vector{Int64}}, total_time::Float64, string_of_ops::Vector; dt=0.01, full=true, enumeration_style::Symbol=:rowwise)
    rows, cols = rows_cols
    N = rows * cols
    couplings = string_of_ops[1:2:end]
    ops = string_of_ops[2:2:end]
    term_ops = [split(op, ",") for op in ops]
    len_op = [length(split(op, ",")) for op in ops]
    trotter_step = round(Int, total_time/dt)
    all_ops = Vector{QuantumOps}()

    for term in 1:length(ops)
        terms = term_ops[term]
        j = len_op[term] - 1
        angle = round(-2 * dt * couplings[term], sigdigits=6)

        if j == 0  # Single-qubit terms
            for row in 1:rows, col in 1:cols
                qubit = _get_qubit_index(rows, cols, row, col, enumeration_style)
                _apply_term(all_ops, terms, angle, [qubit])
            end
        else  # Two-qubit terms
            # Horizontal connections
            for row in 1:rows, col in 1:cols-1
                qubits = [_get_qubit_index(rows, cols, row, c, enumeration_style) for c in col:col+1]
                _apply_term(all_ops, terms, angle, qubits)
            end
            # Vertical connections
            for row in 1:rows-1, col in 1:cols
                qubits = [_get_qubit_index(rows, cols, r, col, enumeration_style) for r in row:row+1]
                _apply_term(all_ops, terms, angle, qubits)
            end
        end
    end

    return full ? vcat(fill(all_ops, trotter_step)...) : all_ops
end