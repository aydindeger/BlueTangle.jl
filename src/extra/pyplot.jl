# # Troubleshooting Plotting Issues

# If you encounter problems with plotting, please follow these steps:

# 1. **Install PyPlot in Julia**: Add the `PyPlot` package to your Julia environment. This package provides an interface to the `matplotlib` library in Python. You can install it using the Julia package manager:
#    ```julia
#    import Pkg
#    Pkg.add("PyPlot")
#    ```

# 2. **Install Python Matplotlib**: Ensure that `matplotlib` is installed in your Python environment. This is a prerequisite for `PyPlot` as it relies on Python's `matplotlib` for plotting. You can install `matplotlib` using `pip`:
#    ```bash
#    pip3 install matplotlib
#    ```

# For detailed documentation and additional information, refer to the [`PyPlot` GitHub page](https://github.com/JuliaPy/plt.jl).

using PyPlot

const AbstractVectorS = Union{AbstractVector,sa.SparseVector}
const AbstractMatrixS = Union{AbstractMatrix,sa.SparseMatrixCSC}

"""
`plot(sample_probs::Vector; rep::Symbol=:int, basis::String="Z")`

Plots the outcome probabilities for quantum measurements.

- `sample_probs`: Vector of tuples, each containing outcome probabilities and number of qubits.

Creates a bar plot showing the probabilities of the most likely outcomes in the specified measurement basis.
"""
function plotq(mVector::Vector{Measurement}, labels::Vector{String}=[""])

    if length(mVector) > 5
        throw("You can only plot five measurements.")
    end

    base_fig_width, base_fig_height = (7, 5)
    x_range = maximum([length(m.bitstr) for m = mVector])
    scale_factor = 0.5
    fig_width = base_fig_width + (x_range * scale_factor)
    fig_width = max(min(fig_width, 16), 5)

    fig, ax = subplots(figsize=(fig_width, base_fig_height), dpi=100)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]

    for (i, r) = enumerate(mVector)

        x_outcome = r.bitstr
        y_sample = r.sample

        wid = 0.8 / i

        lab = labels == [""] ? "($(i-1)) $(r.circuit_name)" : labels[i]
        bars = ax.bar(x_outcome, y_sample, wid, alpha=0.7, color=colors[i], label=lab)

        # Add the probabilities on top of each bar
        for bar in bars
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.005, string(round(height, digits=3)),
                ha="center", va="bottom", color=colors[i])
        end

        # ax.set_xticks(0:x_range)
    end

    # if rep==:bstr
    # ax.set_xlabel("Binary outcomes = left to right (first to last qubit)")
    # else
    ax.set_xlabel("Outcomes")
    # end

    ax.set_ylabel("Probabilities")
    ax.set_title("Outcome Probabilities")
    ax.grid(true)
    ax.legend()
    display(fig)

end



"""
`plot(m::Measurement; rep::Symbol=:int)`

Plots the outcome probabilities for a single quantum measurement.

- `m`: A Measurement object.

Creates a bar plot showing the probabilities of the most likely outcomes from the measurement.
"""
plotq(m::Measurement) = plotq([m])

plotq(state::AbstractVectorS) = plotq([measure(state)])
plotq(states::Vector{<:sa.SparseVector}) = plotq(measure.([states...]))

## circuit drawing

"""
`_draw_gate(ax, op::QuantumOps, pos, gate_width, qubit_lines)`

Internal function to draw a quantum gate on a circuit diagram.

- `ax`: The plot axis.
- `op`: The QuantumOps object representing the quantum operation.
- `pos`: Position on the circuit diagram.
- `gate_width`: Width of the gate drawing.

Draws the specified quantum gate on the given plot axis.
"""
function _draw_gate(ax, op::QuantumOps, pos, gate_width, qubit_lines)

    if isa(op, OpQC)
        c = "tab:red"
        c_t = "tab:red"
    else
        c = "black"
        c_t = "black"
    end

    if isa(op, OpF)
        # Draw a full vertical line for OpF operations
        ax.plot([pos, pos], [-0.2, qubit_lines - 0.8], "black", linewidth=2, markersize=gate_width)
        ax.text(pos, qubit_lines - 0.6, op.name, color="black", ha="center", va="center")
        return
    end

    qubit = op.qubit
    target_qubit = BlueTangle._target_find(op)
    control = BlueTangle._control_find(op)

    if target_qubit > 0
        marker2 = isa(op, OpQC) ? "o" : (BlueTangle._check_swap_invariant(op.mat) ? "o" : "x")
    end

    if control != -2
        # Draw a line for the control-target connection
        ax.plot([pos, pos], [qubit - 1, control - 1], c, markersize=gate_width)
        ax.plot(pos, control - 1, "o", color=c, markersize=gate_width)

        # Draw the control dot
        if op.q == 1
            ax.plot(pos, qubit - 1, "x", color=c, markersize=gate_width)
            if control > qubit
                ax.text(pos, qubit - 1.4, "c-" * op.name, color=c_t, ha="center")
            else
                ax.text(pos, qubit - 0.7, "c-" * op.name, color=c_t, ha="center")
            end
        elseif op.q == 2
            ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width)
            ax.plot(pos, target_qubit - 1, "x", color=c, markersize=gate_width)
            ax.plot([pos, pos], [qubit - 1, target_qubit - 1], c, markersize=gate_width)
            ax.text(pos, target_qubit - 0.8, "c-" * op.name, color=c_t, ha="center")
        end

        # Single qubit gate
    elseif op.q == 1
        ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width)
        ax.text(pos, qubit - 0.7, op.name, color=c_t, ha="center", va="center")

    elseif op.q == 2
        ax.plot([pos, pos], [qubit - 1, target_qubit - 1], c, markersize=gate_width)
        ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width)
        ax.plot(pos, target_qubit - 1, marker2, color=c, markersize=gate_width)
        if target_qubit < qubit
            ax.text(pos, target_qubit - 1.4, op.name, color=c_t, ha="center")
        else
            ax.text(pos, target_qubit - 0.7, op.name, color=c_t, ha="center")
        end
    end
end



"""
    plotq(ops::Vector{QuantumOps}; labels::Vector{String} = [""], fold::Int = 25, reverse_bits::Bool = false)

Plots a quantum circuit diagram from a vector of quantum operations.

# Arguments
- `layers`: Vector of quantum operations or layers of operations
- `labels`: Custom labels for qubits (default: "Qubit (i)")
- `fold`: Number of layers/ops per row before folding to a new row (default: 25, set to 0 to disable folding)
- `reverse_bits`: If true, reverse qubit ordering so last qubit is at top (default: false)

Creates a visual representation of the quantum circuit based on the specified operations and initial qubit states.
When the circuit is longer than `fold` layers, it wraps to multiple rows like Qiskit's matplotlib drawer.
"""
function plotq(layers::Vector; labels::Vector{String}=[""], fold::Int=25, reverse_bits::Bool=false)

    ops = vcat(layers...)

    qubit_lines = maximum([max(op.qubit, BlueTangle._target_find(op), BlueTangle._control_find(op)) for op in ops if !isa(op, OpF)])

    num_ops = length(ops)
    num_layers = length(layers)
    println("layers=$(num_layers), ops=$(num_ops)")

    num = isa(layers, Vector{<:QuantumOps}) ? num_ops : num_layers

    # Calculate folding parameters
    if fold > 0 && num > fold
        num_folds = ceil(Int, num / fold)
        ops_per_fold = fold
    else
        num_folds = 1
        ops_per_fold = num
        fold = num  # No folding needed
    end

    # Drawing constants
    gate_width = 8

    # Calculate figure size with folding
    # Each fold gets its own row of qubit lines - need extra spacing to avoid overlap
    xsize = min(ops_per_fold * 0.6 + 2, 20)  # Width based on ops per fold, max 20
    row_height = qubit_lines * 1.0 + 1.5  # Height per fold row with extra padding
    ysize = num_folds * row_height + 1.0  # Total height

    fig, ax = subplots(figsize=(xsize, ysize), dpi=100)
    ax.axis("off")  # Turn off the axis

    # Draw each fold (row)
    for fold_idx in 1:num_folds
        # Calculate y offset for this fold (rows go from top to bottom)
        y_offset = -(fold_idx - 1) * row_height

        # Calculate the range of operations for this fold
        start_idx = (fold_idx - 1) * fold + 1
        end_idx = min(fold_idx * fold, num)
        ops_in_fold = end_idx - start_idx + 1

        # Draw horizontal lines for qubits in this fold
        for i in 1:qubit_lines
            # Reverse qubit display order if requested
            display_q = reverse_bits ? (qubit_lines - i + 1) : i
            y_pos = (i - 1) + y_offset
            ax.hlines(y_pos, -0.5, ops_in_fold - 0.5, colors="black", linewidth=1)

            # Only show labels on first fold, or on all folds if circuit is long
            label_text = labels == [""] ? "q$display_q" : "$(labels[display_q])"
            ax.text(-0.7, y_pos, label_text, ha="right", va="center", fontsize=9)
        end

        # Add fold indicator on the right if not the last fold
        if fold_idx < num_folds
            x_end = ops_in_fold - 0.3
            for i in 1:qubit_lines
                y_pos = (i - 1) + y_offset
                ax.plot(x_end + 0.1, y_pos, ">", color="gray", markersize=6)
            end
            ax.text(x_end + 0.3, y_offset + qubit_lines / 2 - 0.5, "↓", fontsize=12, color="gray", ha="center", va="center")
        end

        # Draw gates for this fold
        if isa(layers, Vector{<:QuantumOps})
            for global_pos in start_idx:end_idx
                local_pos = global_pos - start_idx  # Position within this fold
                op = ops[global_pos]
                _draw_gate_folded(ax, op, local_pos, gate_width, qubit_lines, y_offset, reverse_bits)
            end
        else
            for global_layer_idx in start_idx:end_idx
                local_pos = global_layer_idx - start_idx  # Position within this fold
                layer = layers[global_layer_idx]
                for op in layer
                    _draw_gate_folded(ax, op, local_pos, gate_width, qubit_lines, y_offset, reverse_bits)
                end
            end
        end
    end

    # Set overall plot limits
    ax.set_xlim(-1.5, ops_per_fold + 0.5)
    ax.set_ylim(-num_folds * row_height, qubit_lines - 0.5)

    tight_layout()
    display(fig)
end

"""
Internal function to draw a gate with y_offset for folded circuits.
"""
function _draw_gate_folded(ax, op::QuantumOps, pos, gate_width, qubit_lines, y_offset, reverse_bits::Bool=false)

    # Helper to map qubit index for display
    map_q(q) = reverse_bits ? (qubit_lines - q + 1) : q

    if isa(op, OpQC)
        c = "tab:red"
        c_t = "tab:red"
    else
        c = "black"
        c_t = "black"
    end

    if isa(op, OpF)
        # Draw a full vertical line for OpF operations
        ax.plot([pos, pos], [y_offset - 0.2, y_offset + qubit_lines - 0.8], "black", linewidth=2, markersize=gate_width)
        ax.text(pos, y_offset + qubit_lines - 0.6, op.name, color="black", ha="center", va="center", fontsize=8)
        return
    end

    qubit = map_q(op.qubit)
    target_qubit_raw = BlueTangle._target_find(op)
    control_raw = BlueTangle._control_find(op)
    target_qubit = target_qubit_raw > 0 ? map_q(target_qubit_raw) : target_qubit_raw
    control = control_raw > 0 ? map_q(control_raw) : control_raw

    if target_qubit_raw > 0
        marker2 = isa(op, OpQC) ? "o" : (BlueTangle._check_swap_invariant(op.mat) ? "o" : "x")
    end

    if control != -2
        # Draw a line for the control-target connection
        ax.plot([pos, pos], [qubit - 1 + y_offset, control - 1 + y_offset], c, markersize=gate_width)
        ax.plot(pos, control - 1 + y_offset, "o", color=c, markersize=gate_width)

        # Draw the control dot
        if op.q == 1
            ax.plot(pos, qubit - 1 + y_offset, "x", color=c, markersize=gate_width)
            label_y = control > qubit ? qubit - 1.3 + y_offset : qubit - 0.6 + y_offset
            ax.text(pos, label_y, "c" * op.name, color=c_t, ha="center", fontsize=7)
        elseif op.q == 2
            ax.plot(pos, qubit - 1 + y_offset, "o", color=c, markersize=gate_width)
            ax.plot(pos, target_qubit - 1 + y_offset, "x", color=c, markersize=gate_width)
            ax.plot([pos, pos], [qubit - 1 + y_offset, target_qubit - 1 + y_offset], c, markersize=gate_width)
            ax.text(pos, target_qubit - 0.7 + y_offset, "c" * op.name, color=c_t, ha="center", fontsize=7)
        end

        # Single qubit gate
    elseif op.q == 1
        ax.plot(pos, qubit - 1 + y_offset, "o", color=c, markersize=gate_width)
        ax.text(pos, qubit - 0.6 + y_offset, op.name, color=c_t, ha="center", va="center", fontsize=8)

    elseif op.q == 2
        ax.plot([pos, pos], [qubit - 1 + y_offset, target_qubit - 1 + y_offset], c, markersize=gate_width)
        ax.plot(pos, qubit - 1 + y_offset, "o", color=c, markersize=gate_width)
        ax.plot(pos, target_qubit - 1 + y_offset, marker2, color=c, markersize=gate_width)
        label_y = target_qubit < qubit ? target_qubit - 1.3 + y_offset : target_qubit - 0.6 + y_offset
        ax.text(pos, label_y, op.name, color=c_t, ha="center", fontsize=8)
    end
end


plotq(circuit::Circuit; labels::Vector{String}=[""]) = plotq(circuit.layers; labels=labels)

plotq(ansatz::AnsatzOptions; labels::Vector{String}=[""]) = plotq(ansatz.ops; labels=labels)


"""
`savefigure(name::String)`
saves figure
"""
savefigure(name::String) = savefig(name)


# ============================================
# Text-based Circuit Drawing
# ============================================

"""
    plotq_text(ops::Vector{<:QuantumOps}; fold::Int=80, labels::Vector{String}=[""], reverse_bits::Bool=false)

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
plotq_text(ops)
```
"""
function plotq_text(layers::Vector; fold::Int=80, labels::Vector{String}=[""], reverse_bits::Bool=false)

    ops = vcat(layers...)

    # Determine if input is flat ops or layered
    is_layered = !isa(layers, Vector{<:QuantumOps})

    # Get number of qubits
    qubit_lines = 0
    for op in ops
        if !isa(op, OpF)
            qubit_lines = max(qubit_lines, op.qubit)
            tq = BlueTangle._target_find(op)
            ctrl = BlueTangle._control_find(op)
            tq > 0 && (qubit_lines = max(qubit_lines, tq))
            ctrl > 0 && (qubit_lines = max(qubit_lines, ctrl))
        end
    end

    num = is_layered ? length(layers) : length(ops)

    # Gate width in characters (for gate name + borders)
    gate_char_width = 5

    # Build the circuit as a 2D character grid
    # Each column represents a layer/op position
    # Each row represents a qubit (with extra rows for labels and spacing)

    # Characters
    wire = '─'
    vline = '│'
    cross = '┼'
    ctrl_dot = '■'
    target_x = '╳'
    box_tl = '┌'
    box_tr = '┐'
    box_bl = '└'
    box_br = '┘'
    box_l = '│'
    box_r = '│'

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
        # 3 rows per qubit: top border, content, bottom (for spacing)
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
    target_qubit_raw = BlueTangle._target_find(op)
    target_qubit = target_qubit_raw > 0 ? map_q(target_qubit_raw) : target_qubit_raw
    control_raw = BlueTangle._control_find(op)
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
plotq_text(op::Op; kwargs...) = plotq_text([op]; kwargs...)
plotq_text(circuit::Circuit; kwargs...) = plotq_text(circuit.layers; kwargs...)
plotq_text(ansatz::AnsatzOptions; kwargs...) = plotq_text(ansatz.ops; kwargs...)
