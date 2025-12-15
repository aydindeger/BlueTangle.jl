using PyCall
ENV["PYCALL_JL_RUNTIME_PYTHON"] = joinpath(@__DIR__, "BlueTangle.jl", "src", "extra", ".venv", "bin", "python")

# ============================================
# Qiskit Circuit Drawing
# ============================================

"""
    plotqiskit(ops::Vector{<:QuantumOps}, output::String="text"; n::Int=0 filename::String="")

Convert BlueTangle ops to a Qiskit circuit and draw it.

# Arguments
- `ops`: Vector of BlueTangle Op objects
- `n`: Number of qubits (auto-detected if 0)
- `output`: Drawing output format: "text", "mpl" (matplotlib), or "latex"
- `filename`: If provided, save the figure to this file (only for "mpl" output)

# Returns
- For "text": returns the ASCII circuit drawing as a string
- For "mpl": displays the matplotlib figure (and saves if filename provided)
- For "latex": returns the LaTeX source

# Example
# In your Julia session after including src.jl
ops = [Op("H", 1), Op("CX", 1, 2), Op("T", 2), Op("CZ", 2, 3)]
plotqiskit(ops)  # Returns text drawing
# For matplotlib figure
plotqiskit(ops; output="mpl")
# Save to file
plotqiskit(ops; output="mpl", filename="circuit.png")
# Works with OpF too
plotqiskit(my_opf)
"""
function plotqiskit(ops::Vector{<:QuantumOps}, output::String="text", filename::String=""; n::Int=0)

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
            if name == "CX" || name == "CNOT"
                qc.ccx(ctrl, q, tq)
            elseif name == "CZ"
                qc.ccz(ctrl, q, tq)
            elseif name == "X"
                qc.ccx(ctrl, q, q)  # This would be CX with extra control
            elseif name == "Z"
                qc.ccz(ctrl, q, q)
            else
                # For other controlled gates, try adding control
                @warn "Controlled gate $name may not be supported directly"
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
            elseif name == "CSX"
                qc.csx(q, tq)
            elseif occursin("CRX", name)
                # Extract angle from CRX(angle)
                m = match(r"CRX\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.crx(angle, q, tq)
                end
            elseif occursin("CRY", name)
                m = match(r"CRY\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.cry(angle, q, tq)
                end
            elseif occursin("CRZ", name)
                m = match(r"CRZ\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.crz(angle, q, tq)
                end
            elseif occursin("RXX", name)
                m = match(r"RXX\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.rxx(angle, q, tq)
                end
            elseif occursin("RYY", name)
                m = match(r"RYY\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.ryy(angle, q, tq)
                end
            elseif occursin("RZZ", name)
                m = match(r"RZZ\(([-\d.]+)\)", name)
                if m !== nothing
                    angle = parse(Float64, m.captures[1])
                    qc.rzz(angle, q, tq)
                end
            else
                @warn "Two-qubit gate $name not directly supported, skipping"
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
            m = match(r"RX\(([-\d.]+)\)", name)
            if m !== nothing
                angle = parse(Float64, m.captures[1])
                qc.rx(angle, q)
            end
        elseif occursin("RY", name)
            m = match(r"RY\(([-\d.]+)\)", name)
            if m !== nothing
                angle = parse(Float64, m.captures[1])
                qc.ry(angle, q)
            end
        elseif occursin("RZ", name)
            m = match(r"RZ\(([-\d.]+)\)", name)
            if m !== nothing
                angle = parse(Float64, m.captures[1])
                qc.rz(angle, q)
            end
        elseif occursin("P", name) || occursin("PHASE", name)
            m = match(r"P(?:HASE)?\(([-\d.]+)\)", name)
            if m !== nothing
                angle = parse(Float64, m.captures[1])
                qc.p(angle, q)
            end
        elseif occursin("U1", name)
            m = match(r"U1\(([-\d.]+)\)", name)
            if m !== nothing
                angle = parse(Float64, m.captures[1])
                qc.p(angle, q)  # U1 = P gate in modern Qiskit
            end
        elseif occursin("U2", name)
            m = match(r"U2\(([-\d.]+),([-\d.]+)\)", name)
            if m !== nothing
                phi = parse(Float64, m.captures[1])
                lam = parse(Float64, m.captures[2])
                qc.u(π / 2, phi, lam, q)
            end
        elseif occursin("U3", name) || occursin("U(", name)
            m = match(r"U[3]?\(([-\d.]+),([-\d.]+),([-\d.]+)\)", name)
            if m !== nothing
                theta = parse(Float64, m.captures[1])
                phi = parse(Float64, m.captures[2])
                lam = parse(Float64, m.captures[3])
                qc.u(theta, phi, lam, q)
            end
        else
            @warn "Gate $name not recognized, skipping"
        end
    end

    # Draw the circuit
    if output == "text"
        return qc.draw(output="text")
    elseif output == "mpl"
        fig = qc.draw(output="mpl")
        if !isempty(filename)
            fig.savefig(filename, dpi=150, bbox_inches="tight")
            println("Saved to $filename")
        end
        return fig
    elseif output == "latex"
        return qc.draw(output="latex_source")
    else
        return qc.draw(output=output)
    end
end

# Convenience method for single Op
plotqiskit(op::Op; kwargs...) = plotqiskit([op]; kwargs...)

# Method for OpF that contains vector of ops
function plotqiskit(opf::OpF; kwargs...)
    if isa(opf.data, Vector)
        return plotqiskit(opf.data; kwargs...)
    else
        error("OpF does not contain a vector of ops")
    end
end