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

# For detailed documentation and additional information, refer to the [`PyPlot` GitHub page](https://github.com/JuliaPy/PyPlot.jl).


"""
`plot_measurement(mVector::Vector{Measurement}; rep::Symbol=:int)`

Plots the outcome probabilities of quantum measurements from a vector of Measurement objects.

- `mVector`: Vector of Measurement objects.
- `rep`: (Optional) Representation of outcomes. Default is `:int`.

Creates a bar plot showing the probabilities of the most likely outcomes.
"""
function plot_measurement(mVector::Vector{Measurement};rep::Symbol=:int)

    sample_prob=[]
    basis_list=[]
    for m = mVector
        push!(sample_prob,(m.int_basis,m.sample))
    end

    plot_measurement(sample_prob,mVector[1].number_of_qubits;rep=rep)

end

"""
`plot_measurement(sample_prob::Vector; rep::Symbol=:int, basis::String="Z")`

Plots the outcome probabilities for quantum measurements.

- `sample_prob`: Vector of tuples, each containing outcome probabilities and number of qubits.
- `rep`: (Optional) Representation of outcomes. Default is `:int`.
- `basis`: (Optional) Measurement basis. Default is "Z".

Creates a bar plot showing the probabilities of the most likely outcomes in the specified measurement basis.
"""
function plot_measurement(sample_prob::Vector,N::Int;rep::Symbol=:int)

    fig, ax = subplots()

    colors=["red","blue","green","orange","black"]

    for (i,r) = enumerate(sample_prob)
    
        r=(bit_to(r[1],N,rep),r[2])

        width=.2i#rand()/2
        bars=ax.bar(r..., width, alpha=0.6, label="$(i)", color=colors[i])

            # Add the probabilities on top of each bar
        for bar in bars
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.005, string(round(height, digits=3)), 
                    ha="center", va="bottom", color=colors[i])
        end

        ax.set_xticks(0:maximum(r[1]))

    end

    if rep==:bstr
        ax.set_xlabel("Binary outcomes = left to right (first to last qubit)")
    else
        ax.set_xlabel("Outcomes")
    end

    ax.set_ylabel("Probabilities")
    ax.set_title("Outcome Probabilities")
    ax.grid(true)
    ax.legend()
    display(fig)

end

"""
`plot_measurement(m::Measurement; rep::Symbol=:int)`

Plots the outcome probabilities for a single quantum measurement.

- `m`: A Measurement object.
- `rep`: (Optional) Representation of outcomes. Default is `:int`.

Creates a bar plot showing the probabilities of the most likely outcomes from the measurement.
"""
plot_measurement(m::Measurement;rep::Symbol=:int)=plot_measurement((m.int_basis,m.sample),m.number_of_qubits;rep=rep,basis=m.measurement_basis)

"""
`plot_measurement(sample_prob::Tuple, N::Int; rep::Symbol=:int, basis::String="Z")`

Plots the outcome probabilities for quantum measurements based on a tuple of sample probabilities and number of qubits.

- `sample_prob`: A tuple containing outcome probabilities and corresponding states.
- `N`: Number of qubits.
- `rep`: (Optional) Representation of outcomes. Default is `:int`.
- `basis`: (Optional) Measurement basis. Default is "Z".

Creates a bar plot showing the probabilities of the most likely outcomes.
"""
function plot_measurement(sample_prob::Tuple,N::Int;rep::Symbol=:int,basis::String="Z")

    fig,ax=subplots(1,1,figsize=(7,4),dpi=100)

    sample_prob=(bit_to(sample_prob[1],N,rep),sample_prob[2])

    bars = ax.bar(sample_prob..., color="red", alpha=0.5)

    # Add the probabilities on top of each bar
    for bar in bars
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height, string(round(height, digits=3)), 
                ha="center", va="bottom")
    end

    ax.set_xticks(0:maximum(sample_prob[1]))

    if rep==:bstr
        ax.set_xlabel("Binary outcomes = left to right (first to last qubit)")
    else
        ax.set_xlabel("Outcomes")
    end

    ax.set_ylabel("Probabilities")
    ax.set_title("Outcome Probabilities at $(basis) basis")
    ax.grid(true)
    # ax.legend()
    display(fig)
end

"""
`_draw_gate(ax, op::QuantumOps, pos, gate_width)`

Internal function to draw a quantum gate on a circuit diagram.

- `ax`: The plot axis.
- `op`: The QuantumOps object representing the quantum operation.
- `pos`: Position on the circuit diagram.
- `gate_width`: Width of the gate drawing.

Draws the specified quantum gate on the given plot axis.
"""
function _draw_gate(ax, op::QuantumOps, pos, gate_width)
    qubit = op.qubit
    
    target_qubit = _target_find(op)

    c="black"
    if op.noise!=false && !_is_it_measurement(op.name)
        c_t="red"
    else
        c_t="black"
    end

    # Single qubit gate
    if op.q == 1
        # Draw the gate symbol (e.g., a circle)

        if _is_it_measurement(op.name)
            ax.plot(pos, qubit - 1, marker=">", color=c, markersize=gate_width*30)
        else
            ax.plot(pos, qubit - 1, "o", color=c, markersize=gate_width*20)
        end
        # Add the gate name
        ax.text(pos, qubit - .8, op.name, color=c_t, ha="center", va="center")

    elseif op.q == 2
        # Draw a line for the control-target connection
        ax.plot([pos, pos], [qubit - 1, target_qubit - 1], c)

        # Draw the control dot
        ax.plot(pos, qubit - 1, "s", color=c)

        if op.name=="CX"
            ax.plot(pos, target_qubit - 1, "x", color=c, markersize=gate_width*20)
        else
            ax.plot(pos, target_qubit - 1, "s", color=c, markersize=gate_width*20)
        end

        if target_qubit<qubit
            ax.text(pos, target_qubit - 1.3, op.name, color=c_t, ha="center")
        else
            ax.text(pos, target_qubit - 0.8, op.name, color=c_t, ha="center")
        end
    end
end

"""
`plot_circuit(circuit::Circuit)`

Plots a quantum circuit diagram from a Circuit object.

- `circuit`: A Circuit object representing the quantum circuit.

Creates a visual representation of the quantum circuit.
"""
plot_circuit(circuit::Circuit)=plot_circuit(circuit.ops)


_target_find(op::QuantumOps)=typeof(op)==ifOp ? -1 : op.target_qubit

"""
`plot_circuit(ops::Vector{QuantumOps}; list_of_initial_qubits::Vector{Int}=Int[])`

Plots a quantum circuit diagram from a vector of quantum operations.

- `ops`: Vector of QuantumOps representing the operations in the circuit.
- `list_of_initial_qubits`: (Optional) List of initial qubit states.

Creates a visual representation of the quantum circuit based on the specified operations and initial qubit states.
"""
function plot_circuit(ops::Vector{QuantumOps};list_of_initial_qubits::Vector{Int}=Int[])

    if isempty(list_of_initial_qubits)
        qubit_lines = maximum([max(op.qubit, _target_find(op)) for op in ops])
    else
        qubit_lines = length(list_of_initial_qubits)
    end

    num_ops = length(ops)

    xsize=num_ops<16 ? 10 : .7num_ops
    ysize=num_ops>5 ? 7 : 1.1qubit_lines

    # Create a new figure and axis
    
    fig, ax = subplots(figsize=(xsize,ysize),dpi=150)
    ax.axis("off")  # Turn off the axis

    # # Set the background color of the axes
    # color=(0.29411764705882354, 0.35294117647058826, 0.5803921568627451, .2)
    # # ax.set_facecolor(color)  # Use any color you like here

    # # Optionally, set the background color of the figure
    # fig.patch.set_facecolor(color)  # Use any color you like here
    
    # Drawing constants
    gate_width = 0.3

    # Set plot limits
    ax.set_ylim(-1, qubit_lines)
    ax.set_xlim(-0.5, num_ops)

    # Draw the horizontal lines for the qubits and label them
    for i in 1:qubit_lines
        ax.hlines(i-1, -0.5, num_ops-0.5, colors="black")
        # Label each qubit line
        if list_of_initial_qubits ==Int[]
            ax.text(-0.7, i-1, "Qubit $i", ha="right", va="center")
        else
            ax.text(-0.7, i-1, "Qubit $(i) [$(list_of_initial_qubits[i])]", ha="right", va="center")
        end
    end

    # Add gates to the circuit
    for (pos, op) in enumerate(ops)
        _draw_gate(ax, op, pos-1, gate_width)  # Position gates with some offset
    end

    # Show the plot
    display(fig)

    # return fig
end