using PyPlot
import BlueTangle: plotq, savefigure

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

    ax.set_xticks(sort(unique(vcat([m.bitstr for m in mVector]...))))
    ax.set_ylim(0, 1.15 * maximum(maximum(m.sample) for m in mVector))
    ax.set_ylabel("Probabilities")
    ax.set_title("Outcome Probabilities")
    ax.grid(true)
    ax.legend()
    display(fig)
    return fig

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
function plotq(layers::Union{Vector{<:QuantumOps},Vector{<:AbstractVector{<:QuantumOps}}}; labels::Vector{String}=[""], fold::Int=25, reverse_bits::Bool=false)

    ops = QuantumOps[op for group in layers for op in (group isa QuantumOps ? (group,) : group)]

    qubit_lines = get_N(ops)

    num_ops = length(ops)
    num_layers = length(layers)

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
    ax.set_ylim(-(num_folds - 1) * row_height - 0.5, qubit_lines - 0.5)

    tight_layout()
    display(fig)
    return fig
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


plotq(circuit::Circuit; kwargs...) = plotq(circuit.layers; kwargs...)

plotq(ansatz::AnsatzOptions; kwargs...) = plotq(ansatz.ops; kwargs...)


"""
`savefigure(name::String)`
saves figure
"""
savefigure(name::String) = savefig(name)


# Keep the historical helper name while sharing the package's text renderer.
plotq_text(args...; kwargs...) = plotqt(args...; kwargs...)
