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

"""
`plot(sample_probs::Vector; rep::Symbol=:int, basis::String="Z")`

Plots the outcome probabilities for quantum measurements.

- `sample_probs`: Vector of tuples, each containing outcome probabilities and number of qubits.

Creates a bar plot showing the probabilities of the most likely outcomes in the specified measurement basis.
"""
function plotq(mVector::Vector{Measurement},labels::Vector{String}=[""])

    if length(mVector)>5
        throw("You can only plot five measurements.")
    end

    base_fig_width, base_fig_height = (7, 5)
    x_range = maximum([length(m.bitstr) for m=mVector])
    scale_factor = 0.5
    fig_width = base_fig_width + (x_range * scale_factor)
    fig_width = max(min(fig_width, 16), 5)

    fig, ax = subplots(figsize=(fig_width, base_fig_height),dpi=100)

    colors=["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown"]

    for (i,r) = enumerate(mVector)

        x_outcome=r.bitstr
        y_sample=r.sample

        wid=.8/i

        lab=labels==[""] ? "($(i-1)) $(r.circuit_name)" : labels[i]
        bars=ax.bar(x_outcome, y_sample, wid, alpha=0.7, color=colors[i], label=lab)

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
plotq(m::Measurement)=plotq([m])

plotq(state::sa.SparseVector)=plotq([measure(state)])
plotq(states::Vector{<:sa.SparseVector})=plotq(measure.([states...]))

## circuit drawing

_target_find(op::QuantumOps)=typeof(op)==ifOp ? -1 : op.target_qubit

function _control_find(op::QuantumOps)
    if isa(op,ifOp) || isa(op,QC)
        return -2
    else
        return op.control
    end
end

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

    if isa(op, QC)
        c = "tab:red"
        c_t = "tab:red"
    else
        c = "black"
        c_t = "black"
    end

    if isa(op, OpF)
        # Draw a full vertical line for OpF operations
        ax.plot([pos, pos], [-.2, qubit_lines-.8], "black", linewidth=2)
        ax.text(pos, qubit_lines - 0.6, op.name, color="black", ha="center", va="center")
        return
    end

    qubit = op.qubit
    target_qubit = _target_find(op)
    control = _control_find(op)

    if target_qubit>0
        marker2=isa(op,QC) ? "o" : (BlueTangle._swap_control_target(op.mat)==op.mat ? "x" : "o")
    end

    if control != -2
        # Draw a line for the control-target connection
        ax.plot([pos, pos], [qubit - 1, control - 1], c)
        ax.plot(pos, control - 1, "o", color=c, markersize=gate_width * 20)

        # Draw the control dot
        if op.q == 1
            ax.plot(pos, qubit - 1, "x", color=c, markersize=gate_width * 20)
            if control > qubit
                ax.text(pos, qubit - 1.4, "c-" * op.name, color=c_t, ha="center")
            else
                ax.text(pos, qubit - 0.7, "c-" * op.name, color=c_t, ha="center")
            end
        elseif op.q == 2
            ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width * 20)
            ax.plot(pos, target_qubit - 1, "x", color=c, markersize=gate_width * 20)
            ax.plot([pos, pos], [qubit - 1, target_qubit - 1], c)
            ax.text(pos, target_qubit - 0.8, "c-" * op.name, color=c_t, ha="center")
        end

    # Single qubit gate
    elseif op.q == 1
        ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width * 20)
        ax.text(pos, qubit - 0.7, op.name, color=c_t, ha="center", va="center")

    elseif op.q == 2
        ax.plot([pos, pos], [qubit - 1, target_qubit - 1], c, markersize=gate_width * 20)
        ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width * 20)
        ax.plot(pos, target_qubit - 1, marker2, color=c, markersize=gate_width * 20)
        if target_qubit < qubit
            ax.text(pos, target_qubit - 1.4, op.name, color=c_t, ha="center")
        else
            ax.text(pos, target_qubit - 0.7, op.name, color=c_t, ha="center")
        end
    end
end



"""
    plot(ops::Vector{QuantumOps}; labels::Vector{String} = [""])

Plots a quantum circuit diagram from a vector of quantum operations.

Creates a visual representation of the quantum circuit based on the specified operations and initial qubit states.
"""
function plotq(layers::Vector; labels::Vector{String} = [""])
    ops = vcat(layers...)

    qubit_lines = maximum([max(op.qubit, _target_find(op), _control_find(op)) for op in ops if !isa(op, OpF)])

    num_ops = length(ops)
    num_layers = length(layers)
    println("layers=$(num_layers), ops=$(num_ops)")

    num = isa(layers, Vector{<:QuantumOps}) ? num_ops : num_layers

    # Adjust the xsize calculation
    xsize = num * 1.5  # Provide horizontal space based on number of operations
    ysize = qubit_lines * 1  # Vertical space based on qubit lines
    
    fig, ax = subplots(figsize=(xsize, ysize), dpi=100)
    ax.axis("off")  # Turn off the axis
    
    # Drawing constants
    gate_width = 0.4

    # Set plot limits
    ax.set_ylim(-1, qubit_lines)
    # ax.set_xlim(-0.5, num_ops)

    # Draw the horizontal lines for qubits and label them
    for i in 1:qubit_lines
        ax.hlines(i - 1, -0.5, num - 0.5, colors="black")  # Horizontal line
        label_text = labels == [""] ? "Qubit ($i)" : "$(labels[i]) ($i)"
        ax.text(-0.7, i - 1, label_text, ha="right", va="center")
    end

    if isa(layers, Vector{<:QuantumOps})
        for (pos, op) in enumerate(ops)
            _draw_gate(ax, op, pos - 1, gate_width, qubit_lines)  # Position gates with some offset
        end
    else
        for (layer_idx, layer) in enumerate(layers)
            for op in layer
                _draw_gate(ax, op, layer_idx - 1, gate_width, qubit_lines)  # Use layer_idx as horizontal position
            end
        end
    end

    display(fig)
end


plotq(circuit::Circuit; labels::Vector{String} = [""])=plotq(circuit.layers;labels=labels)

plotq(ansatz::AnsatzOptions; labels::Vector{String} = [""])=plotq(ansatz.ops;labels=labels)


"""
`savefigure(name::String)`
saves figure
"""
savefigure(name::String)=savefig(name)
