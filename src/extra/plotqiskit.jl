using PyCall
ENV["PYCALL_JL_RUNTIME_PYTHON"] = joinpath(@__DIR__, "BlueTangle.jl", "src", "extra", ".venv", "bin", "python")


# ============================================
# Qiskit Circuit Drawing
# ============================================

"""
    plotqiskit(ops::Vector{<:QuantumOps}, output::String="mpl"; n::Int=0, fold=nothing, filename::String="")

Convert BlueTangle ops to a Qiskit circuit and draw it.

# Arguments
- `ops`: Vector of BlueTangle Op objects
- `n`: Number of qubits (auto-detected if 0)
- `output`: Drawing output format: "text", "mpl" (matplotlib), or "latex"
- `fold`: Optional Qiskit fold value. If omitted, uses 1.5x Qiskit's default fold (25 -> 38).
- `filename`: If provided, save the figure to this file (only for "mpl" output)

# Returns
- For "text": returns the ASCII circuit drawing as a string
- For "mpl": displays the matplotlib figure (and saves if filename provided)
- For "latex": returns the LaTeX source

# Example
# In your Julia session after including src.jl
ops = [Op("H", 1), Op("CX", 1, 2), Op("T", 2), Op("CZ", 2, 3)]
plotqiskit(ops)  # Returns matplotlib figure (default)
# For text drawing
plotqiskit(ops, "text")
# Save to file
plotqiskit(ops, "mpl", "circuit.png")
# Works with OpF too
plotqiskit(my_opf)
"""
function plotqiskit(
    ops::Vector{<:QuantumOps},
    output::String="mpl",
    filename::String="";
    n::Int=0,
    fold::Union{Nothing,Int}=nothing,
)

    # Auto-detect number of qubits if not provided
    if n == 0
        for op in ops
            if isa(op, Op)
                n = max(n, op.qubit)
                op.target_qubit > 0 && (n = max(n, op.target_qubit))
                op.control > 0 && (n = max(n, op.control))
            end
        end
    end

    # Create Qiskit circuit via PyCall
    qiskit = pyimport("qiskit")
    QuantumCircuit = qiskit.QuantumCircuit

    qc = QuantumCircuit(n)

    _parse_angle_expr(expr::AbstractString) = begin
        ex = replace(strip(expr), "π" => "pi", "PI" => "pi")
        try
            val = eval(Meta.parse(ex))
            val isa Number || error("not numeric")
            Float64(val)
        catch
            error("Failed to parse gate angle expression '$expr' in plotqiskit.")
        end
    end

    # Gate name mapping from BlueTangle to Qiskit
    # BlueTangle uses 1-indexed qubits, Qiskit uses 0-indexed
    for op in ops
        if !isa(op, Op)
            continue  # Skip non-Op quantum operations (OpF, ifOp, etc.)
        end

        name = uppercase(op.name)
        q = op.qubit - 1  # Convert to 0-indexed
        tq = op.target_qubit > 0 ? op.target_qubit - 1 : -1
        ctrl = op.control > 0 ? op.control - 1 : -1

        # Handle controlled gates (CCX, CCZ, etc.)
        if ctrl >= 0
            if name == "CX" || name == "CNOT" || name == "CCX"
                qc.ccx(ctrl, q, tq)
            elseif name == "CZ" || name == "CCZ"
                qc.ccz(ctrl, q, tq)
            elseif name == "X"
                qc.ccx(ctrl, q, q)  # This would be CX with extra control
            elseif name == "Z"
                qc.ccz(ctrl, q, q)
            else
                error("Unsupported controlled gate '$name' in plotqiskit.")
            end
            continue
        end

        # Two-qubit gates
        if tq >= 0
            if name == "CX" || name == "CNOT"
                qc.cx(q, tq)
            elseif name == "CZ"
                qc.cz(q, tq)
            elseif name == "CY"
                qc.cy(q, tq)
            elseif name == "SWAP"
                qc.swap(q, tq)
            elseif name == "ISWAP"
                qc.iswap(q, tq)
            elseif name == "CH"
                qc.ch(q, tq)
            elseif name == "CS"
                qc.cs(q, tq)
            elseif name == "CSD" || name == "CSDG" || name == "CSDAG"
                if hasproperty(qc, :csdg)
                    qc.csdg(q, tq)
                else
                    qc.cp(-π / 2, q, tq)
                end
            elseif name == "CT"
                qc.cp(π / 4, q, tq)
            elseif name == "CTD" || name == "CTDG" || name == "CTDAG"
                qc.cp(-π / 4, q, tq)
            elseif name == "CSX"
                qc.csx(q, tq)
            elseif occursin("CRX", name)
                # Extract angle from CRX(angle)
                m = match(r"CRX\(([^)]+)\)", name)
                m === nothing && error("Could not parse CRX parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.crx(angle, q, tq)
            elseif occursin("CRY", name)
                m = match(r"CRY\(([^)]+)\)", name)
                m === nothing && error("Could not parse CRY parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.cry(angle, q, tq)
            elseif occursin("CRZ", name)
                m = match(r"CRZ\(([^)]+)\)", name)
                m === nothing && error("Could not parse CRZ parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.crz(angle, q, tq)
            elseif occursin("RXX", name)
                m = match(r"RXX\(([^)]+)\)", name)
                m === nothing && error("Could not parse RXX parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.rxx(angle, q, tq)
            elseif occursin("RYY", name)
                m = match(r"RYY\(([^)]+)\)", name)
                m === nothing && error("Could not parse RYY parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.ryy(angle, q, tq)
            elseif occursin("RZZ", name)
                m = match(r"RZZ\(([^)]+)\)", name)
                m === nothing && error("Could not parse RZZ parameters from '$name' in plotqiskit.")
                angle = _parse_angle_expr(m.captures[1])
                qc.rzz(angle, q, tq)
            else
                error("Unsupported two-qubit gate '$name' in plotqiskit.")
            end
            continue
        end

        # Single-qubit gates
        if name == "X"
            qc.x(q)
        elseif name == "Y"
            qc.y(q)
        elseif name == "Z"
            qc.z(q)
        elseif name == "H"
            qc.h(q)
        elseif name == "S"
            qc.s(q)
        elseif name == "SD" || name == "SDG" || name == "SDAG"
            qc.sdg(q)
        elseif name == "T"
            qc.t(q)
        elseif name == "TD" || name == "TDG" || name == "TDAG"
            qc.tdg(q)
        elseif name == "SX"
            qc.sx(q)
        elseif name == "SXDG"
            qc.sxdg(q)
        elseif name == "I" || name == "ID"
            qc.id(q)
        elseif occursin("RX", name)
            m = match(r"RX\(([^)]+)\)", name)
            m === nothing && error("Could not parse RX parameters from '$name' in plotqiskit.")
            angle = _parse_angle_expr(m.captures[1])
            qc.rx(angle, q)
        elseif occursin("RY", name)
            m = match(r"RY\(([^)]+)\)", name)
            m === nothing && error("Could not parse RY parameters from '$name' in plotqiskit.")
            angle = _parse_angle_expr(m.captures[1])
            qc.ry(angle, q)
        elseif occursin("RZ", name)
            m = match(r"RZ\(([^)]+)\)", name)
            m === nothing && error("Could not parse RZ parameters from '$name' in plotqiskit.")
            angle = _parse_angle_expr(m.captures[1])
            qc.rz(angle, q)
        elseif occursin("P", name) || occursin("PHASE", name)
            m = match(r"P(?:HASE)?\(([^)]+)\)", name)
            m === nothing && error("Could not parse P/PHASE parameters from '$name' in plotqiskit.")
            angle = _parse_angle_expr(m.captures[1])
            qc.p(angle, q)
        elseif occursin("U1", name)
            m = match(r"U1\(([^)]+)\)", name)
            m === nothing && error("Could not parse U1 parameters from '$name' in plotqiskit.")
            angle = _parse_angle_expr(m.captures[1])
            qc.p(angle, q)  # U1 = P gate in modern Qiskit
        elseif occursin("U2", name)
            m = match(r"U2\(([^,]+),([^)]+)\)", name)
            m === nothing && error("Could not parse U2 parameters from '$name' in plotqiskit.")
            phi = _parse_angle_expr(m.captures[1])
            lam = _parse_angle_expr(m.captures[2])
            qc.u(π / 2, phi, lam, q)
        elseif occursin("U3", name) || occursin("U(", name)
            m = match(r"U[3]?\(([^,]+),([^,]+),([^)]+)\)", name)
            m === nothing && error("Could not parse U/U3 parameters from '$name' in plotqiskit.")
            theta = _parse_angle_expr(m.captures[1])
            phi = _parse_angle_expr(m.captures[2])
            lam = _parse_angle_expr(m.captures[3])
            qc.u(theta, phi, lam, q)
        else
            error("Unsupported single-qubit gate '$name' in plotqiskit.")
        end
    end

    # Qiskit default fold is 25 for mpl; use 1.5x when user doesn't specify.
    fold_use = fold === nothing ? 50 : fold

    # Draw the circuit
    if output == "text"
        return qc.draw(output="text", fold=fold_use)
    elseif output == "mpl"
        fig = qc.draw(output="mpl", fold=fold_use)
        if !isempty(filename)
            fig.savefig(filename, dpi=150, bbox_inches="tight")
            println("Saved to $filename")
        end
        return fig
    elseif output == "latex"
        return qc.draw(output="latex_source", fold=fold_use)
    else
        return qc.draw(output=output, fold=fold_use)
    end
end

# Convenience method for single Op
plotqiskit(op::Op; kwargs...) = plotqiskit([op]; kwargs...)

# Convenience method for QASM file path or QASM string
function plotqiskit(
    qasm_input::AbstractString,
    output::String="mpl",
    filename::String="";
    n::Int=0,
    fold::Union{Nothing,Int}=nothing,
)
    ops = if isfile(qasm_input)
        from_qasm_file(qasm_input)
    elseif occursin("OPENQASM", uppercase(qasm_input))
        from_qasm(qasm_input)
    else
        error("String input to plotqiskit must be a QASM file path or OPENQASM string.")
    end
    return plotqiskit(ops, output, filename; n=n, fold=fold)
end

# Method for OpF that contains vector of ops
function plotqiskit(opf::OpF; kwargs...)
    if isa(opf.data, Vector)
        return plotqiskit(opf.data; kwargs...)
    else
        error("OpF does not contain a vector of ops")
    end
end


"""
    from_qasm(qasm_content::AbstractString) -> Vector{Op}

Import OpenQASM text and convert it to BlueTangle `Op` objects.
"""
function from_qasm(qasm_content::AbstractString)
    qc = try
        qasm2 = pyimport("qiskit.qasm2")
        pycall(qasm2.loads, PyObject, qasm_content)
    catch
        qiskit = pyimport("qiskit")
        qiskit.QuantumCircuit.from_qasm_str(qasm_content)
    end

    gate_finder = _load_gate_finder_bridge(force_reload=false)
    raw = gate_finder._qiskit_circuit_to_gate_tuples(qc)
    normalized = [gate_finder._normalize_gate_tuple(g) for g in raw]
    return _py_encoding_circuit_to_ops(normalized)
end

from_qasm_file(path::AbstractString) = from_qasm(read(path, String))