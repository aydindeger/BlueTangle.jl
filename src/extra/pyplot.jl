# Loaded by the PyPlot package extension. The core text renderer supplies the
# circuit layout, so drawing never evaluates gates or changes the random stream.
using PyPlot
import BlueTangle: _plotq_pyplot, savefigure

const _plot_colors = ["#3267a8", "#d18032", "#448b78", "#9b669b", "#b85e60", "#768947"]

"""
    plotq(measurements, labels=[""]; backend=:pyplot, rep=:bits, max_outcomes=32, show=true)

Draw grouped measurement probabilities. Binary labels put q1 on the left;
`rep=:int` uses integer outcomes. When there are more than `max_outcomes`, retain
the most probable outcomes and combine the remainder in an explicit Other bar.
Signed estimates from measurement mitigation remain visible below zero.
"""
function _plotq_pyplot(measurements::Vector{Measurement}, labels::Vector{String}=[""];
               rep::Symbol=:bits, max_outcomes::Int=32, show::Bool=true)
    tick_labels, values, series_labels, N = BlueTangle._measurement_plot_data(measurements, labels; rep=rep, max_outcomes=max_outcomes)
    positions = collect(0:length(tick_labels)-1)
    fig, ax = subplots(figsize=(clamp(0.55length(tick_labels) + 1.8, 5.4, 16), 3.8), dpi=120)
    fig.set_facecolor("white")
    ax.set_axisbelow(true)
    ax.grid(axis="y", color="#e5e9ef", linewidth=0.8)
    for side in ("top", "right", "left")
        ax.spines[side].set_visible(false)
    end
    ax.spines["bottom"].set_color("#b9c3ce")
    ax.tick_params(axis="both", length=0, labelsize=9, pad=7)
    width = 0.78 / length(measurements)
    lower, upper = 0.0, 0.0
    for (i, measurement) in enumerate(measurements)
        probabilities = values[:, i]
        lower = min(lower, minimum(probabilities; init=0.0))
        upper = max(upper, maximum(probabilities; init=0.0))
        label = series_labels[i]
        offset = (i - (length(measurements) + 1) / 2) * width
        bars = ax.bar(positions .+ offset, probabilities, width=width,
                      color=_plot_colors[mod1(i, length(_plot_colors))], label=label, zorder=3)
        if (length(measurements) == 1 && length(tick_labels) <= 16) ||
           (length(measurements) == 2 && length(tick_labels) <= 4)
            for (bar, probability) in zip(bars, probabilities)
                abs(probability) < 1e-12 && continue
                ax.annotate("$(round(100probability; digits=1))%",
                            (bar.get_x() + bar.get_width()/2, probability),
                            xytext=(0, probability >= 0 ? 4 : -4), textcoords="offset points",
                            ha="center", va=probability >= 0 ? "bottom" : "top", fontsize=8, color="#384657")
            end
        end
    end
    span = max(upper - lower, 0.1)
    ax.set_ylim(lower < 0 ? lower - 0.16span : 0, upper + 0.18span)
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, rotation=(length(tick_labels) > 12 || N > 6) ? 45 : 0,
                       ha=(length(tick_labels) > 12 || N > 6) ? "right" : "center")
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.PercentFormatter(xmax=1))
    signed = any(<(0), values)
    ax.set_ylabel(signed ? "Quasiprobability" : "Probability", fontsize=10, labelpad=10)
    basis = measurements[1].measurement_basis == "0" ? "Z" : measurements[1].measurement_basis
    ax.set_xlabel(rep == :int ? "Outcome ($basis basis)" : "Outcome ($basis basis; q1 at left)", fontsize=10, labelpad=11)
    ax.set_title(signed ? "Measurement estimates" : "Measurement probabilities", loc="left", fontsize=12, fontweight="medium", pad=length(measurements) > 1 ? 36 : 15)
    if length(measurements) > 1
        ax.legend(frameon=false, ncol=min(length(measurements), 3), loc="lower left",
                  bbox_to_anchor=(0, 1.01), borderaxespad=0, fontsize=9)
    end
    fig.tight_layout(pad=1.2)
    show && display(fig)
    return fig
end

_plotq_pyplot(measurement::Measurement; kwargs...) = _plotq_pyplot([measurement]; kwargs...)
_plotq_pyplot(state::Union{AbstractVector{<:Number},AbstractMatrix}; kwargs...) = _plotq_pyplot(measure(state); kwargs...)
_plotq_pyplot(states::Vector{<:AbstractVector{<:Number}}; kwargs...) = _plotq_pyplot(measure.(states); kwargs...)

function _gate_box(ax, x, ys, width, label; face="#eef3fa", edge="#3267a8")
    bottom, top = extrema(ys)
    padding = occursin('\n', label) ? 0.40 : 0.28
    rectangle = PyPlot.matplotlib.patches.FancyBboxPatch(
        (x-width/2, bottom-padding), width, top-bottom+2padding,
        boxstyle="round,pad=0.025,rounding_size=0.06", facecolor=face,
        edgecolor=edge, linewidth=1.15, zorder=3)
    ax.add_patch(rectangle)
    ax.text(x, (bottom+top)/2, label, ha="center", va="center", fontsize=9,
            color="#243347", zorder=4, parse_math=false)
    if length(ys) > 1
        for (port, y) in enumerate(ys)
            ax.text(x-width/2+0.07, y, string(port), ha="left", va="center",
                    fontsize=7, color=edge, zorder=4)
        end
    end
end

_circuit_gate_width(op) = op isa ifOp ?
    max(1.6, 0.085maximum(length, split(BlueTangle._circuit_gate_label(op), ';'))) :
    max(0.65, 0.085length(BlueTangle._circuit_gate_label(op)) + 0.22)

function _draw_circuit_gate(ax, op, x, width, N, yoffset, reverse_bits)
    y(q) = (reverse_bits ? q-1 : N-q) + yoffset
    ink = "#33465a"
    function connect(qubits; xpos=x)
        ys = y.(qubits)
        ax.plot([xpos, xpos], [minimum(ys), maximum(ys)], color=ink, linewidth=1.15, zorder=2)
    end
    dot(q; xpos=x) = ax.plot(xpos, y(q), "o", color=ink, markersize=5, zorder=5)
    function target(q)
        # Marker sizes use points, keeping ⊕ circular at any figure aspect ratio.
        ax.plot(x, y(q), marker="o", markersize=13, markerfacecolor="white",
                markeredgecolor=ink, markeredgewidth=1.3, linestyle="none", zorder=3)
        ax.plot(x, y(q), marker="+", markersize=11, color=ink,
                markeredgewidth=1.1, linestyle="none", zorder=4)
    end
    symbol = BlueTangle._circuit_symbol(op)
    if symbol in (:cx, :cz, :swap)
        qubits = [op.qubit, op.target_qubit]
        op.control > 0 && push!(qubits, op.control)
        connect(qubits)
        if symbol == :swap
            for q in (op.qubit, op.target_qubit)
                ax.plot(x, y(q), "x", color=ink, markersize=8, markeredgewidth=1.5, zorder=4)
            end
        else
            dot(op.qubit)
            symbol == :cx ? target(op.target_qubit) : dot(op.target_qubit)
        end
        op.control > 0 && dot(op.control)
        return
    end
    if op isa Op && op.q == 1 && op.control > 0 && uppercase(op.name) == "X" && op.mat == gate.X
        connect([op.qubit, op.control]); dot(op.control); target(op.qubit)
        return
    end

    qubits = if op isa Op
        op.q == 1 ? [op.qubit] : [op.qubit, op.target_qubit]
    else
        collect(BlueTangle._operated_qubits(op))
    end
    isempty(qubits) && (qubits = collect(1:N))
    label = BlueTangle._circuit_gate_label(op)
    if op isa ifOp
        label = "$(op.name) on q$(op.qubit)\n0: $(BlueTangle._circuit_branch_label(op.if01[1]))\n1: $(BlueTangle._circuit_branch_label(op.if01[2]))"
    end
    face, edge = op isa OpQC ? ("#fcf0e2", "#b77b37") :
                 op isa ifOp ? ("#e9f4ef", "#438b74") :
                 op isa OpF ? ("#f1edf7", "#8a6ca6") : ("#eef3fa", "#3267a8")
    _gate_box(ax, x, y.(qubits), width, label; face=face, edge=edge)
    if op isa Op && op.control > 0
        connector = length(qubits) == 1 ? x : x-width/2-0.14
        connect([qubits; op.control]; xpos=connector)
        dot(op.control; xpos=connector)
    end
end

"""
    plotq(ops_or_layers; backend=:pyplot, labels=[""], fold=25, reverse_bits=false, nqubits=nothing, show=true)

Draw a circuit using gate boxes and conventional controlled-X, controlled-Z and
SWAP symbols. Qubit 1 is at the top. `fold` counts visual columns; zero disables
folding. Returns the figure for further customization or explicit saving.
"""
function _plotq_pyplot(input::AbstractVector;
               labels::Vector{String}=[""], fold::Int=25, reverse_bits::Bool=false,
               nqubits::Union{Nothing,Int}=nothing, show::Bool=true)
    fold >= 0 || throw(ArgumentError("fold must be nonnegative"))
    columns, N = BlueTangle._circuit_plot_data(input; nqubits=nqubits)
    names = BlueTangle._circuit_qubit_labels(labels, N)
    if isempty(columns)
        fig, ax = subplots(figsize=(4, 1.3), dpi=120)
        ax.axis("off")
        ax.text(0.5, 0.5, "Empty circuit", ha="center", va="center", color="#657487")
        show && display(fig)
        return fig
    end
    widths = [maximum(_circuit_gate_width, column) for column in columns]
    # Balance folded rows instead of leaving a nearly empty final row.
    row_count = fold == 0 ? 1 : cld(length(columns), fold)
    step = cld(length(columns), row_count)
    ranges = [first:min(first+step-1, length(columns)) for first in 1:step:length(columns)]
    row_widths = [sum(widths[range]) + 0.36length(range) for range in ranges]
    row_height = max(N, 1) + 0.85
    fig, ax = subplots(figsize=(max(4, 0.75maximum(row_widths) + 1.0), max(1.7, 0.52row_height * length(ranges))), dpi=120)
    fig.set_facecolor("white")
    ax.axis("off")
    for (row, range) in enumerate(ranges)
        offset = -(row-1)*row_height
        ending = row_widths[row]
        for q in 1:N
            ypos = (reverse_bits ? q-1 : N-q) + offset
            ax.plot([-0.23, ending], [ypos, ypos], color="#aab7c5", linewidth=0.9, zorder=1)
            ax.text(-0.42, ypos, names[q], ha="right", va="center", fontsize=10, color="#344459", parse_math=false)
        end
        x = 0.0
        for index in range
            width = widths[index]
            center = x + width/2 + 0.12
            for op in columns[index]
                _draw_circuit_gate(ax, op, center, _circuit_gate_width(op), N, offset, reverse_bits)
            end
            x += width + 0.36
        end
        if row < length(ranges)
            ax.text(ending+0.04, offset+(N-1)/2, "›", ha="left", va="center", fontsize=16, color="#7e8c9c")
        end
    end
    ax.set_xlim(-0.85 - 0.04maximum(length, names), maximum(row_widths)+0.25)
    ax.set_ylim(-(length(ranges)-1)*row_height-0.55, N-0.45)
    fig.tight_layout(pad=0.5)
    show && display(fig)
    return fig
end

_plotq_pyplot(op::QuantumOps; kwargs...) = _plotq_pyplot([op]; kwargs...)
_plotq_pyplot(ansatz::AnsatzOptions; kwargs...) = _plotq_pyplot(ansatz.ops; kwargs...)

"""
    savefigure(path; kwargs...)
    savefigure(figure, path; kwargs...)

Save the current or specified figure, trimming unused outer margins by default.
"""
function savefigure(path::AbstractString; bbox_inches="tight", kwargs...)
    PyPlot.savefig(path; bbox_inches=bbox_inches, kwargs...)
    return path
end
function savefigure(figure::PyPlot.Figure, path::AbstractString; bbox_inches="tight", kwargs...)
    figure.savefig(path; bbox_inches=bbox_inches, kwargs...)
    return path
end

plotq_text(args...; kwargs...) = plotqt(args...; kwargs...)
