# ============================================
# Text-based Circuit Drawing
# ============================================

"""
    plotqt(ops::Vector{<:QuantumOps}; fold::Int=80, labels::Vector{String}=[""], reverse_bits::Bool=false)

Generate a text-based circuit diagram similar to Qiskit's text output.

# Arguments
- `ops`: Vector of quantum operations or layers
- `fold`: Maximum width before folding to next row (default: 80 characters)
- `labels`: Custom labels for qubits
- `reverse_bits`: If true, reverse qubit ordering so last qubit is at top (default: false)

# Returns
- Prints the ASCII circuit diagram

# Example
```julia
ops = [Op("H", 1), Op("CX", 1, 2), Op("T", 2)]
plotqt(ops)
```
"""
function plotqt(layers::Vector; fold::Int=80, labels::Vector{String}=[""], reverse_bits::Bool=false)

    ops = vcat(layers...)

    # Determine if input is flat ops or layered
    is_layered = !isa(layers, Vector{<:QuantumOps})

    # Get number of qubits
    qubit_lines = 0
    for op in ops
        if !isa(op, OpF)
            qubit_lines = max(qubit_lines, op.qubit)
            tq = _target_find(op)
            ctrl = _control_find(op)
            tq > 0 && (qubit_lines = max(qubit_lines, tq))
            ctrl > 0 && (qubit_lines = max(qubit_lines, ctrl))
        end
    end

    num = is_layered ? length(layers) : length(ops)

    # Gate width in characters (for gate name + borders)
    gate_char_width = 5

    # Characters
    wire = '─'

    # Calculate label width
    if labels == [""]
        max_label_len = length("q_$(qubit_lines): ")
    else
        max_label_len = maximum(length.(labels)) + 3
    end

    # Each layer takes gate_char_width characters
    total_width = num * gate_char_width

    # Check if folding is needed
    usable_width = fold - max_label_len - 3  # Account for labels and fold marker
    ops_per_fold = max(1, usable_width ÷ gate_char_width)

    if num > ops_per_fold && fold > 0
        num_folds = ceil(Int, num / ops_per_fold)
    else
        num_folds = 1
        ops_per_fold = num
    end

    result = IOBuffer()

    for fold_idx in 1:num_folds
        start_idx = (fold_idx - 1) * ops_per_fold + 1
        end_idx = min(fold_idx * ops_per_fold, num)
        fold_width = (end_idx - start_idx + 1) * gate_char_width

        # Initialize grid for this fold
        rows_per_qubit = 2
        grid_height = qubit_lines * rows_per_qubit
        grid_width = fold_width

        grid = fill(' ', grid_height, grid_width)

        # Fill in wire lines (every other row starting at row 1)
        for q in 1:qubit_lines
            row = (q - 1) * rows_per_qubit + 1
            for col in 1:grid_width
                grid[row, col] = wire
            end
        end

        # Process operations in this fold
        if is_layered
            for layer_idx in start_idx:end_idx
                local_pos = layer_idx - start_idx
                col_start = local_pos * gate_char_width + 1
                layer = layers[layer_idx]
                for op in layer
                    _draw_text_gate!(grid, op, col_start, gate_char_width, rows_per_qubit, qubit_lines, reverse_bits)
                end
            end
        else
            for op_idx in start_idx:end_idx
                local_pos = op_idx - start_idx
                col_start = local_pos * gate_char_width + 1
                op = ops[op_idx]
                _draw_text_gate!(grid, op, col_start, gate_char_width, rows_per_qubit, qubit_lines, reverse_bits)
            end
        end

        # Print this fold
        for q in 1:qubit_lines
            # Map qubit for display (reverse if requested)
            display_q = reverse_bits ? (qubit_lines - q + 1) : q
            # Print label
            if labels == [""]
                label = "q_$(display_q-1): "  # 0-indexed like Qiskit
            else
                label = "$(labels[display_q]): "
            end
            print(result, lpad(label, max_label_len))

            # Print grid rows for this qubit
            row = (q - 1) * rows_per_qubit + 1
            for col in 1:grid_width
                print(result, grid[row, col])
            end

            # Add fold continuation marker if not last fold
            if fold_idx < num_folds
                print(result, "─»")
            end
            println(result)

            # Print spacing row (vertical connections)
            spacing_row = (q - 1) * rows_per_qubit + 2
            if q < qubit_lines
                print(result, " "^max_label_len)
                for col in 1:grid_width
                    print(result, grid[spacing_row, col])
                end
                println(result)
            end
        end

        # Add separator between folds (with extra blank lines for clarity)
        if fold_idx < num_folds
            println(result)
            println(result)  # Extra blank line for spacing
            println(result, "─"^40)  # Separator line
            println(result)
        end
    end

    # Print the result directly (for Jupyter compatibility)
    circuit_str = String(take!(result))
    print(circuit_str)
    return nothing
end

"""
Internal function to draw a gate on the text grid.
"""
function _draw_text_gate!(grid, op::QuantumOps, col_start::Int, gate_width::Int, rows_per_qubit::Int, qubit_lines::Int, reverse_bits::Bool=false)

    # Helper to map qubit index for display
    map_q(q) = reverse_bits ? (qubit_lines - q + 1) : q

    if isa(op, OpF)
        # Draw vertical barrier line
        for q in 1:qubit_lines
            row = (q - 1) * rows_per_qubit + 1
            mid = col_start + gate_width ÷ 2
            if mid <= size(grid, 2)
                grid[row, mid] = '║'
            end
        end
        return
    end

    qubit = map_q(op.qubit)
    target_qubit_raw = _target_find(op)
    target_qubit = target_qubit_raw > 0 ? map_q(target_qubit_raw) : target_qubit_raw
    control_raw = _control_find(op)
    control = control_raw > 0 ? map_q(control_raw) : control_raw

    # Get short gate name (max 3 chars for box)
    name = uppercase(op.name)
    if length(name) > 3
        name = name[1:3]
    end

    row_q = (qubit - 1) * rows_per_qubit + 1
    mid = col_start + gate_width ÷ 2

    if control != -2
        # Controlled gate
        row_ctrl = (control - 1) * rows_per_qubit + 1

        # Draw control dot
        if mid <= size(grid, 2)
            grid[row_ctrl, mid] = '■'
        end

        # Draw vertical connection
        min_row = min(row_q, row_ctrl)
        max_row = max(row_q, row_ctrl)
        for r in min_row:max_row
            if r != row_q && r != row_ctrl && mid <= size(grid, 2)
                if grid[r, mid] == '─'
                    grid[r, mid] = '┼'
                else
                    grid[r, mid] = '│'
                end
            end
        end

        # Draw target
        if op.q == 1
            # Single qubit controlled gate (like CX)
            if mid <= size(grid, 2)
                grid[row_q, mid] = '⊕'
            end
        else
            # Two-qubit controlled gate
            row_tq = (target_qubit - 1) * rows_per_qubit + 1
            if mid <= size(grid, 2)
                grid[row_q, mid] = '■'
                grid[row_tq, mid] = '⊕'
            end
        end

    elseif target_qubit > 0
        # Two-qubit gate (like CX, CZ, SWAP)
        row_tq = (target_qubit - 1) * rows_per_qubit + 1

        min_row = min(row_q, row_tq)
        max_row = max(row_q, row_tq)

        # Draw vertical line between qubits
        for r in min_row:max_row
            if mid <= size(grid, 2)
                if r == row_q || r == row_tq
                    grid[r, mid] = '■'
                elseif grid[r, mid] == '─'
                    grid[r, mid] = '┼'
                else
                    grid[r, mid] = '│'
                end
            end
        end

        # Also fill spacing rows
        for r in (min_row+1):(max_row-1)
            if mid <= size(grid, 2) && r <= size(grid, 1)
                grid[r, mid] = '│'
            end
        end

    else
        # Single qubit gate - draw box with name
        if col_start <= size(grid, 2) && col_start + gate_width - 1 <= size(grid, 2)
            # Draw gate as: ┤name├
            grid[row_q, col_start] = '┤'

            # Center the name
            name_start = col_start + 1
            for (i, c) in enumerate(name)
                if name_start + i - 1 <= size(grid, 2)
                    grid[row_q, name_start+i-1] = c
                end
            end

            grid[row_q, col_start+gate_width-1] = '├'
        end
    end
end

# Convenience methods
plotqt(op::Op; kwargs...) = plotqt([op]; kwargs...)
plotqt(circuit::Circuit; kwargs...) = plotqt(circuit.layers; kwargs...)
plotqt(ansatz::AnsatzOptions; kwargs...) = plotqt(ansatz.ops; kwargs...)
