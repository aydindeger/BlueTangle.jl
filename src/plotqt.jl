# Shared circuit layout for the text and optional graphical renderers.
function _circuit_plot_data(input::AbstractVector; nqubits::Union{Nothing,Int}=nothing)
    groups = if all(op -> op isa QuantumOps, input)
        [QuantumOps[op] for op in input]
    elseif all(group -> group isa AbstractVector && all(op -> op isa QuantumOps, group), input)
        [QuantumOps[group...] for group in input]
    else
        throw(ArgumentError("provide quantum operations or layers of quantum operations"))
    end
    ops = QuantumOps[op for group in groups for op in group]
    required = get_N(ops)
    N = isnothing(nqubits) ? required : nqubits
    N >= required && N >= 0 || throw(ArgumentError("nqubits must cover every operated qubit"))
    # Opaque whole-register functions have no inferable qubit count.
    N == 0 && !isempty(ops) && (N = 1)

    # Disjoint gates can still overlap on paper when a connector crosses an
    # intermediate wire. Split only those collisions into adjacent columns.
    columns = Vector{QuantumOps}[]
    for group in groups
        lanes = Vector{QuantumOps}[]
        occupied = Vector{UnitRange{Int}}[]
        for op in group
            qubits = collect(_operated_qubits(op))
            span = isempty(qubits) ? (1:N) : (minimum(qubits):maximum(qubits))
            lane = findfirst(ranges -> all(r -> isempty(intersect(r, span)), ranges), occupied)
            if isnothing(lane)
                push!(lanes, QuantumOps[op])
                push!(occupied, [span])
            else
                push!(lanes[lane], op)
                push!(occupied[lane], span)
            end
        end
        append!(columns, lanes)
    end
    return columns, N
end

function _circuit_operation_name(op)
    name = replace(op.name, '\n' => ' ')
    op isa Op || return name
    parsed = match(r"^([A-Za-z][A-Za-z0-9]*)\(([^()]*)\)$", name)
    isnothing(parsed) && return name
    uppercase(parsed.captures[1]) in gates_with_phase || return name
    # Shorten literal roundoff in standard gate parameters without evaluating
    # expressions or modifying arbitrary custom labels.
    parameters = tryparse.(Float64, strip.(split(parsed.captures[2], ',')))
    all(p -> !isnothing(p) && isfinite(p), parameters) || return name
    return "$(parsed.captures[1])($(join(round.(parameters; sigdigits=6), ",")))"
end

function _circuit_branch_label(ops)
    isempty(ops) && return "I"
    return join((uppercase(op.name) == "I" ? "I" : "$(_circuit_operation_name(op))@q$(op.qubit)" for op in ops), ", ")
end

function _circuit_gate_label(op::QuantumOps)
    if op isa ifOp
        return "$(op.name)? 0: $(_circuit_branch_label(op.if01[1])); 1: $(_circuit_branch_label(op.if01[2]))"
    elseif op isa OpQC
        name = uppercase(op.name) in ("RES", "RESET") ? "Reset" : replace(op.name, '_' => ' ')
        return isempty(op.type) ? name : "$name($(op.type))"
    end
    return _circuit_operation_name(op)
end

function _circuit_qubit_labels(labels, N)
    if labels == [""] || isempty(labels)
        return ["q$q" for q in 1:N]
    end
    length(labels) == N || throw(ArgumentError("provide one label per qubit"))
    return labels
end

function _circuit_symbol(op::QuantumOps)
    op isa Op && op.q == 2 && op.mat isa AbstractMatrix || return :box
    name = uppercase(op.name)
    name in ("CX", "CNOT") && op.mat == gate.CX && return :cx
    name == "CZ" && op.mat == gate.CZ && return :cz
    name == "SWAP" && op.mat == gate.SWAP && return :swap
    return :box
end

function _text_gate!(grid, op, first_column, width, N, reverse_bits)
    row(q) = 2 * (reverse_bits ? N - q : q - 1) + 1
    middle = first_column + width ÷ 2
    function connect(qubits)
        rows = row.(qubits)
        for r in minimum(rows):maximum(rows)
            grid[r, middle] = isodd(r) ? '┼' : '│'
        end
    end
    function box(q, label)
        text = collect("[" * label * "]")
        start = first_column + (width - length(text)) ÷ 2
        grid[row(q), start:start+length(text)-1] = text
    end
    name = uppercase(op.name)
    if _circuit_symbol(op) in (:cx, :cz, :swap)
        qubits = [op.qubit, op.target_qubit]
        op.control > 0 && push!(qubits, op.control)
        connect(qubits)
        grid[row(op.qubit), middle] = name == "SWAP" ? '×' : '●'
        grid[row(op.target_qubit), middle] = name == "SWAP" ? '×' : name == "CZ" ? '●' : '⊕'
        op.control > 0 && (grid[row(op.control), middle] = '●')
    elseif op isa Op && op.q == 1 && op.control > 0
        connect([op.qubit, op.control])
        grid[row(op.control), middle] = '●'
        name == "X" && op.mat == gate.X ? (grid[row(op.qubit), middle] = '⊕') : box(op.qubit, _circuit_gate_label(op))
    else
        qubits = op isa Op ? (op.q == 1 ? [op.qubit] : [op.qubit, op.target_qubit]) : collect(_operated_qubits(op))
        isempty(qubits) && (qubits = collect(1:N))
        connected = op isa Op && op.control > 0 ? [qubits; op.control] : qubits
        length(connected) > 1 && connect(connected)
        op isa Op && op.control > 0 && (grid[row(op.control), middle] = '●')
        first_qubit = op isa OpF ? first(qubits) : op.qubit
        box(first_qubit, _circuit_gate_label(op))
        for (port, q) in enumerate(qubits)
            q == first_qubit || box(q, op isa ifOp ? "if" : string(port))
        end
    end
end

"""
    plotqt(ops_or_layers; fold=80, labels=[""], reverse_bits=false, nqubits=nothing, io=stdout)

Print a compact circuit with q1 at the top, gate labels, and conventional
controlled-X, controlled-Z and SWAP symbols. `fold` is the approximate line
width in characters; use zero for a single line. For an opaque whole-register
`OpF`, supply `nqubits` when its register size cannot be inferred.
Literal numeric gate parameters are displayed with six significant digits.
"""
function plotqt(input::AbstractVector; fold::Int=80, labels::Vector{String}=[""],
                reverse_bits::Bool=false, nqubits::Union{Nothing,Int}=nothing, io::IO=stdout)
    fold >= 0 || throw(ArgumentError("fold must be nonnegative"))
    columns, N = _circuit_plot_data(input; nqubits=nqubits)
    if isempty(columns)
        println(io, "(empty circuit)")
        return nothing
    end
    names = _circuit_qubit_labels(labels, N)
    label_width = maximum(length, names) + 2
    widths = [max(5, maximum(length(_circuit_gate_label(op)) + 4 for op in column)) for column in columns]
    first_column = 1
    while first_column <= length(columns)
        last_column = first_column
        used = widths[first_column]
        while last_column < length(columns) && (fold == 0 || label_width + used + widths[last_column+1] <= fold)
            last_column += 1
            used += widths[last_column]
        end
        grid = fill(' ', 2N - 1, used)
        grid[1:2:end, :] .= '─'
        start = 1
        for index in first_column:last_column
            for op in columns[index]
                _text_gate!(grid, op, start, widths[index], N, reverse_bits)
            end
            start += widths[index]
        end
        for r in axes(grid, 1)
            if isodd(r)
                q = (r + 1) ÷ 2
                reverse_bits && (q = N - q + 1)
                print(io, lpad(names[q] * ": ", label_width))
            else
                print(io, " "^label_width)
            end
            println(io, rstrip(String(grid[r, :])))
        end
        last_column < length(columns) && println(io)
        first_column = last_column + 1
    end
    return nothing
end

plotqt(op::QuantumOps; kwargs...) = plotqt([op]; kwargs...)
plotqt(circuit::Circuit; nqubits::Int=circuit.stats.N, kwargs...) = plotqt(circuit.layers; nqubits=nqubits, kwargs...)
plotqt(ansatz::AnsatzOptions; kwargs...) = plotqt(ansatz.ops; kwargs...)

# Optional graphics extend this function; the public default remains text even
# after PyPlot is loaded.
function _plotq_pyplot end

function _measurement_plot_data(measurements, labels; rep=:bits, max_outcomes=32)
    isempty(measurements) && throw(ArgumentError("provide at least one measurement"))
    rep in (:bits, :bstr, :int) || throw(ArgumentError("rep must be :bits or :int"))
    max_outcomes > 0 || throw(ArgumentError("max_outcomes must be positive"))
    labels == [""] || length(labels) == length(measurements) ||
        throw(ArgumentError("provide one label per measurement"))
    outcomes = sort(unique(vcat([collect(m.bitstr) for m in measurements]...)))
    distributions = [Dict(zip(m.bitstr, m.sample)) for m in measurements]
    omitted = eltype(outcomes)[]
    if length(outcomes) > max_outcomes
        importance = [sum(abs(get(d, outcome, 0.0)) for d in distributions) for outcome in outcomes]
        selected = sortperm(importance; rev=true)[1:max_outcomes]
        omitted = outcomes[setdiff(eachindex(outcomes), selected)]
        outcomes = sort(outcomes[selected])
    end
    N = maximum(m.number_of_qubits for m in measurements)
    ticks = rep == :int ? string.(outcomes) : [string(outcome; base=2, pad=N) for outcome in outcomes]
    isempty(omitted) || push!(ticks, "Other")
    values = zeros(Float64, length(ticks), length(measurements))
    for (j, distribution) in enumerate(distributions)
        for (i, outcome) in enumerate(outcomes)
            values[i, j] = get(distribution, outcome, 0.0)
        end
        isempty(omitted) || (values[end, j] = sum(get(distribution, outcome, 0.0) for outcome in omitted))
    end
    names = labels == [""] ? [m.circuit_name for m in measurements] : labels
    return ticks, values, names, N
end

_plotq_text(args...; kwargs...) = plotqt(args...; kwargs...)
_plotq_text(measurement::Measurement; kwargs...) = _plotq_text([measurement]; kwargs...)
_plotq_text(state::Union{AbstractVector{<:Number},AbstractMatrix}; kwargs...) = _plotq_text(measure(state); kwargs...)
_plotq_text(states::Vector{<:AbstractVector{<:Number}}; kwargs...) = _plotq_text(measure.(states); kwargs...)

function _plotq_text(measurements::Vector{Measurement}, labels::Vector{String}=[""];
                     rep::Symbol=:bits, max_outcomes::Int=32, io::IO=stdout)
    ticks, values, names, _ = _measurement_plot_data(measurements, labels; rep=rep, max_outcomes=max_outcomes)
    basis = measurements[1].measurement_basis == "0" ? "Z" : measurements[1].measurement_basis
    title = any(<(0), values) ? "Measurement estimates (quasiprobabilities)" : "Measurement probabilities"
    println(io, "$title ($basis basis", rep == :int ? ")" : "; q1 at left)")
    outcome_width = max(7, maximum(length, ticks; init=0))
    if length(measurements) == 1
        scale = maximum(abs, values; init=0.0)
        for (i, tick) in enumerate(ticks)
            value = values[i, 1]
            count = scale == 0 ? 0 : round(Int, 24abs(value)/scale)
            bar = (value < 0 ? "−" : " ") * "█"^count
            println(io, rpad(tick, outcome_width+2), rpad(bar, 27), "$(round(100value; digits=2))%")
        end
    else
        widths = [max(11, length(name)+2) for name in names]
        println(io, rpad("Outcome", outcome_width+2), join(rpad(name, width) for (name, width) in zip(names, widths)))
        for (i, tick) in enumerate(ticks)
            probabilities = [rpad("$(round(100values[i,j]; digits=2))%", widths[j]) for j in eachindex(names)]
            println(io, rpad(tick, outcome_width+2), join(probabilities))
        end
    end
    return nothing
end

"""
    plotq(input; backend=:text, kwargs...)

Show a circuit or measurement without extra dependencies. The default text
backend uses gate boxes, conventional controls, and q1 at the top (or the left
of outcome bitstrings). `plotqt` is the direct circuit-text interface.

For a graphical figure, load `PyPlot` and request `backend=:pyplot` explicitly.
Graphical calls return a figure that can be saved with `savefigure(fig, path)`.
"""
function plotq(args...; backend::Symbol=:text, kwargs...)
    backend == :text && return _plotq_text(args...; kwargs...)
    if backend == :pyplot
        isempty(methods(_plotq_pyplot)) &&
            throw(ArgumentError("graphical plotting is optional; install PyPlot and run `using PyPlot` first"))
        return _plotq_pyplot(args...; kwargs...)
    end
    throw(ArgumentError("backend must be :text or :pyplot"))
end
