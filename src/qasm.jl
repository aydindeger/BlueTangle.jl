###############
# OpenQASM I/O
###############

_qi0(i::Int) = i - 1  # 1→0

"""
    to_qasm(ops::Vector{<:QuantumOps};
            nqubits::Union{Int,Nothing}=nothing,
            qreg::AbstractString="q",
            creg::AbstractString="c") -> String

Emit OpenQASM 2.0 for a vector of `Op`s, including parametrized gates
(RX/RY/RZ, U1/U2/U3, CP, RXX/RYY/RZZ/RXY, GIVENS, FSIM, SWAPA) and controls.
"""
function to_qasm(ops::Vector{<:QuantumOps};
                 nqubits::Union{Int,Nothing}=nothing,
                 qreg::AbstractString="q",
                 creg::AbstractString="c")::String
    isempty(ops) && throw(ArgumentError("ops vector is empty"))

    N = isnothing(nqubits) ? get_N(ops) : nqubits
    N < 1 && throw(ArgumentError("nqubits must be ≥ 1"))
    N >= get_N(ops) || throw(ArgumentError("nqubits is smaller than an operation's qubit index."))
    all(name -> occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", name), (qreg, creg)) ||
        throw(ArgumentError("Invalid QASM register name."))
    qreg != creg || throw(ArgumentError("Quantum and classical register names must differ."))

    # collect opaque decls we need to emit
    opaques = Set{String}()

    # op → QASM lines
    function _op_to_qasm(op::QuantumOps)::Vector{String}
        # Special-case: reset is represented internally as OpQC("res", ...) (γ=1)
        if op isa OpQC
            nameU = uppercase(op.name)
            if nameU == "RES" || nameU == "RESET"
                qi = _qi0(op.qubit)
                return ["reset $qreg[$qi];"]
            end
            throw(ArgumentError("Cannot export quantum channel '$(op.name)' to OpenQASM; only RES/RESET (γ=1) is supported."))
        end

        # Everything else: export only plain Op
        op isa Op || throw(ArgumentError("Unsupported operation for export: $(typeof(op))."))
        op = op::Op

        q = op.qubit; t = op.target_qubit; c = op.control
        qi = _qi0(q)
        ti = t >= 1 ? _qi0(t) : t
        ci = c >= 1 ? _qi0(c) : c

        nameU = uppercase(strip(op.name))
        base = _clean_name(op.name)
        args = _name_with_arg_bool(base) ? join(_get_name_and_arguments(op.name)[3], ",") : nothing
        op.mat isa Function && throw(ArgumentError("Bind variational parameters before QASM export."))
        op.mat ≈ gates(op.name) || throw(ArgumentError("Custom matrix for $(op.name) cannot be exported by name."))

        # measurements
        if nameU == "MZ" || nameU == "M(Z)"
            return ["measure $qreg[$qi] -> $creg[$qi];"]
        elseif nameU == "MX" || nameU == "M(X)"
            return ["h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];", "h $qreg[$qi];"]
        elseif nameU == "MY" || nameU == "M(Y)"
            return ["sdg $qreg[$qi];", "h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];", "h $qreg[$qi];", "s $qreg[$qi];"]
        end

        # reset
        if nameU == "RES" || nameU == "RESET"
            return ["reset $qreg[$qi];"]
        end

        # 1-qubit (incl. parametrized) + optional control
        if t < 0
            if c >= 1
                # ---- Controlled 1q gates ----
                if base == "X";  return ["cx $qreg[$ci],$qreg[$qi];"] end
                if base == "Y";  return ["cy $qreg[$ci],$qreg[$qi];"] end
                if base == "Z";  return ["cz $qreg[$ci],$qreg[$qi];"] end
                if base == "H";  return ["ch $qreg[$ci],$qreg[$qi];"] end
                if base == "I";  return ["id $qreg[$qi];"] end
                if base in ("S", "SD", "T", "TD")
                    angle = base == "S" ? pi/2 : base == "SD" ? -pi/2 : base == "T" ? pi/4 : -pi/4
                    return ["cu1($angle) $qreg[$ci],$qreg[$qi];"]
                end

                # Controlled rotations: export as crx/cry/crz (NOT cu1/cu3)
                if base == "RX"
                    push!(opaques, "opaque crx(theta) a,b;")
                    return ["crx($(args)) $qreg[$ci],$qreg[$qi];"]
                elseif base == "RY"
                    push!(opaques, "opaque cry(theta) a,b;")
                    return ["cry($(args)) $qreg[$ci],$qreg[$qi];"]
                elseif base == "RZ"
                    push!(opaques, "opaque crz(theta) a,b;")
                    return ["crz($(args)) $qreg[$ci],$qreg[$qi];"]
                end

                # Phase/U-gates keep cu1/cu3 forms
                if base == "P" || base == "U1"
                    return ["cu1($(args)) $qreg[$ci],$qreg[$qi];"]
                elseif base == "U2"
                    return ["cu3(pi/2,$(args)) $qreg[$ci],$qreg[$qi];"]
                elseif base == "U3"
                    return ["cu3($(args)) $qreg[$ci],$qreg[$qi];"]
                else
                    throw(ArgumentError("Controlled $(op.name) not supported"))
                end
            else
                # plain 1q (no control)
                if base == "P" || base == "U1"; return ["u1($(args)) $qreg[$qi];"] end
                if base == "U2";                return ["u2($(args)) $qreg[$qi];"] end
                if base == "U3";                return ["u3($(args)) $qreg[$qi];"] end
                if base == "RX" || base == "RY" || base == "RZ"
                    return [lowercase(base) * "(" * (args === nothing ? "" : args) * ") $qreg[$qi];"]
                end
                if     nameU == "I";      return ["id $qreg[$qi];"] end
                if     nameU == "X";      return ["x $qreg[$qi];"] end
                if     nameU == "Y";      return ["y $qreg[$qi];"] end
                if     nameU == "Z";      return ["z $qreg[$qi];"] end
                if     nameU == "H";      return ["h $qreg[$qi];"] end
                if     nameU == "S";      return ["s $qreg[$qi];"] end
                if     nameU == "SD";     return ["sdg $qreg[$qi];"] end
                if     nameU == "T";      return ["t $qreg[$qi];"] end
                if     nameU == "TD";     return ["tdg $qreg[$qi];"] end
                if     nameU == "SX" || nameU == "XSQRT"; return ["sx $qreg[$qi];"] end
                if     nameU == "SXDG" || nameU == "SXD" || nameU == "XSQRTDG" || nameU == "XSQRTD"; return ["sxdg $qreg[$qi];"] end
                throw(ArgumentError("Unsupported 1-qubit op: $(op.name)"))
            end
        end

        # Only these two-qubit operators have an exported extra control.
        c >= 1 && !(base in ("CX", "CNOT", "CY", "CZ", "SWAP")) &&
            throw(ArgumentError("Controlled $(op.name) is not supported by QASM export."))
        if base == "CX" || base == "CNOT"
            if c >= 1; return ["ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];"] else return ["cx $qreg[$qi],$qreg[$ti];"] end
        elseif base == "CY"
            if c >= 1
                return ["sdg $qreg[$ti];", "ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];", "s $qreg[$ti];"]
            else
                return ["cy $qreg[$qi],$qreg[$ti];"]
            end
        elseif base == "CZ"
            if c >= 1
                return ["h $qreg[$ti];", "ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];", "h $qreg[$ti];"]
            else
                return ["cz $qreg[$qi],$qreg[$ti];"]
            end
        elseif base == "SWAP"
            if c >= 1; return ["cswap $qreg[$ci],$qreg[$qi],$qreg[$ti];"] else return ["swap $qreg[$qi],$qreg[$ti];"] end
        elseif base == "ISWAP"
            push!(opaques, "opaque iswap a,b;"); return ["iswap $qreg[$qi],$qreg[$ti];"]
        elseif base == "FSWAP"
            push!(opaques, "opaque fswap a,b;"); return ["fswap $qreg[$qi],$qreg[$ti];"]
        elseif base == "SYC"
            push!(opaques, "opaque syc a,b;");    return ["syc $qreg[$qi],$qreg[$ti];"]
        elseif base == "ECR"
            push!(opaques, "opaque ecr a,b;");    return ["ecr $qreg[$qi],$qreg[$ti];"]
        end

        if base == "CP"
            return ["cu1($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "RXX"
            push!(opaques, "opaque rxx(theta) a,b;"); return ["rxx($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "RYY"
            push!(opaques, "opaque ryy(theta) a,b;"); return ["ryy($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "RZZ"
            push!(opaques, "opaque rzz(theta) a,b;"); return ["rzz($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "RXY"
            push!(opaques, "opaque rxy(theta) a,b;"); return ["rxy($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "GIVENS"
            push!(opaques, "opaque givens(theta) a,b;"); return ["givens($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "FSIM"
            push!(opaques, "opaque fsim(theta,phi) a,b;"); return ["fsim($(args)) $qreg[$qi],$qreg[$ti];"]
        elseif base == "SWAPA"
            push!(opaques, "opaque swapa(alpha) a,b;"); return ["swapa($(args)) $qreg[$qi],$qreg[$ti];"]
        end

        throw(ArgumentError("Unsupported operation for export: $(op.name)"))
    end

    header = String["OPENQASM 2.0;", "include \"qelib1.inc\";"]
    decls  = String["qreg $qreg[$N];", "creg $creg[$N];"]
    body   = String[]
    for op in ops
        append!(body, _op_to_qasm(op))
    end
    opaque_lines = sort!(collect(opaques))
    return join(vcat(header, opaque_lines, decls, body), '\n') * '\n'
end

"""
    to_qasm(ops::Vector{<:QuantumOps}, filename::AbstractString;
            nqubits::Union{Int,Nothing}=nothing,
            qreg::AbstractString="q",
            creg::AbstractString="c") -> String

Like `to_qasm(ops; kwargs...)`, but also writes the emitted OpenQASM to `filename`.
"""
function to_qasm(ops::Vector{<:QuantumOps}, filename::AbstractString; kwargs...)::String
    isempty(filename) && throw(ArgumentError("filename must be a non-empty string"))
    qasm = to_qasm(ops; kwargs...)
    open(filename, "w") do io
        write(io, qasm)
    end
    return qasm
end

"""
    to_qasm(circuit::Circuit; kwargs...) -> String
    to_qasm(circuit::Circuit, filename::AbstractString; kwargs...) -> String

Convenience overloads that export a compiled `Circuit` by flattening `circuit.layers`.
"""
to_qasm(circuit::Circuit; kwargs...)::String = to_qasm(vcat(circuit.layers...); kwargs...)
to_qasm(circuit::Circuit, filename::AbstractString; kwargs...)::String =
    to_qasm(vcat(circuit.layers...), filename; kwargs...)


"""
    from_qasm(qasm_or_path::AbstractString; is_file=nothing) -> Vector{QuantumOps}

Read the supported OpenQASM 2 gate subset from text or a file. Quantum registers
are flattened in declaration order, preserving idle qubits. Numeric parameters
use the same restricted arithmetic parser as `gates`. Indexed operands and
whole-register broadcasts are supported. Set `is_file=true` or `false` to select
file or string input explicitly. Parse errors identify the source and line.
Custom gate definitions and classical `if` statements must be expanded before
import. Measurement records follow operation order, not classical-register layout.
"""
function from_qasm(qasm_or_path::AbstractString; is_file::Union{Nothing,Bool}=nothing)::Vector{QuantumOps}
    text = _qasm_to_source(qasm_or_path; is_file=is_file)
    source_name = is_file === true || (is_file !== false && isfile(qasm_or_path)) ? abspath(qasm_or_path) : "<string>"
    source = replace(text, r"//[^\n]*" => "")
    isempty(strip(source)) || endswith(rstrip(source), ';') ||
        throw(ArgumentError("$source_name:$(count(==('\n'),source)+1): Missing ';' after the final statement."))
    statements = Tuple{Int,String}[]
    line = 1
    for chunk in split(source, ';')
        first_line = line + count(==('\n'), match(r"^\s*",chunk).match)
        line += count(==('\n'),chunk)
        isempty(strip(chunk)) || push!(statements,(first_line,String(strip(chunk))))
    end
    registers = Dict{String,UnitRange{Int}}()
    classical_registers = Dict{String,Int}()
    nqubits = 0
    implicit_register = nothing
    ops = QuantumOps[]

    function qubit(token)
        m = match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]$", strip(token))
        m === nothing && throw(ArgumentError("Expected an indexed qubit, got: $token"))
        name, index = m.captures[1], parse(Int, m.captures[2])
        if isempty(registers)
            if implicit_register === nothing
                implicit_register = name
            elseif implicit_register != name
                throw(ArgumentError("Declare quantum registers before using multiple register names."))
            end
            return index + 1
        end
        haskey(registers, name) || throw(ArgumentError("Unknown quantum register: $name"))
        index < length(registers[name]) || throw(ArgumentError("Qubit $token is outside its register."))
        return first(registers[name]) + index
    end

    function operand(token; classical=false)
        name = strip(token)
        if classical
            if haskey(classical_registers,name)
                return collect(1:classical_registers[name])
            end
            m = match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]$",name)
            m === nothing && throw(ArgumentError("Expected a classical register or indexed bit, got: $token"))
            register, index = m.captures[1], parse(Int,m.captures[2])
            if !isempty(classical_registers)
                haskey(classical_registers,register) || throw(ArgumentError("Unknown classical register: $register"))
                index < classical_registers[register] || throw(ArgumentError("Classical bit $token is outside its register."))
            elseif !isempty(registers)
                throw(ArgumentError("Declare the classical destination register before measurement."))
            end
            return [index+1]
        end
        return haskey(registers,name) ? collect(registers[name]) : [qubit(name)]
    end

    for (line,statement) in statements
        try
            occursin(r"^(?:gate|if)\b"i,statement) &&
                throw(ArgumentError("Expand custom gates and classical conditions before importing QASM."))
            if occursin(r"^OPENQASM\s+2(?:\.0)?$"i, statement) ||
               occursin(r"^include\s+\"qelib1\.inc\"$"i, statement) ||
               occursin(r"^(?:opaque|barrier)\s+"i, statement)
                continue
            end
            if (m = match(r"^(qreg|creg)\s+([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]$", statement)) !== nothing
                kind, name, count = m.captures[1], m.captures[2], parse(Int, m.captures[3])
                isempty(ops) || throw(ArgumentError("Declare quantum registers before operations."))
                count > 0 || throw(ArgumentError("Register size must be positive."))
                (haskey(registers,name) || haskey(classical_registers,name)) && throw(ArgumentError("Duplicate register: $name"))
                if kind == "qreg"
                    registers[name] = nqubits + 1:nqubits + count
                    nqubits += count
                else
                    classical_registers[name] = count
                end
                continue
            end
            if (m = match(r"^measure\s+(.+?)\s*->\s*(.+)$"i, statement)) !== nothing
                sources = operand(m.captures[1])
                destinations = operand(m.captures[2]; classical=true)
                length(sources) == length(destinations) || throw(DimensionMismatch("Measurement registers must have the same size."))
                append!(ops,[Op("MZ",q) for q in sources])
                continue
            end
            if (m = match(r"^reset\s+(.+)$"i, statement)) !== nothing
                append!(ops,[Op("RES",q) for q in operand(m.captures[1])])
                continue
            end

            m = match(r"^([A-Za-z][A-Za-z0-9]*)\s*(?:\((.*)\))?\s+(.+)$", statement)
            m === nothing && throw(ArgumentError("Unsupported QASM statement: $statement"))
            name = lowercase(m.captures[1])
            arguments = m.captures[2]
            operands = operand.(split(m.captures[3], ','))
            width = maximum(length,operands)
            all(qs -> length(qs) in (1,width),operands) || throw(DimensionMismatch("Broadcast register sizes must agree."))
            parameterized = arguments !== nothing
            for index in 1:width
                qubits = [qs[length(qs) == 1 ? 1 : index] for qs in operands]
                if !parameterized && name in ("id", "i", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg")
                    length(qubits) == 1 || throw(ArgumentError("$name requires one qubit."))
                    aliases = Dict("id"=>"I", "sdg"=>"SD", "tdg"=>"TD", "sxdg"=>"SXDG")
                    push!(ops, Op(get(aliases, name, uppercase(name)), only(qubits)))
                elseif parameterized && name in ("rx", "ry", "rz", "u", "u1", "u2", "u3", "p")
                    length(qubits) == 1 || throw(ArgumentError("$name requires one qubit."))
                    gate_name = name == "u" ? "U3" : name
                    push!(ops, Op("$gate_name($arguments)", only(qubits)))
                elseif name in ("cx", "cnot", "cy", "cz", "ch", "crx", "cry", "crz", "cu1", "cu2", "cu3", "cp")
                    length(qubits) == 2 || throw(ArgumentError("$name requires two qubits."))
                    has_parameters = name in ("crx", "cry", "crz", "cu1", "cu2", "cu3", "cp")
                    parameterized == has_parameters || throw(ArgumentError("Wrong parameter form for $name."))
                    base = name == "cnot" ? "X" : name == "cp" ? "P" : uppercase(name[2:end])
                    gate_name = parameterized ? "$base($arguments)" : base
                    push!(ops, Op(gate_name, qubits[2]; control=qubits[1]))
                elseif !parameterized && name in ("ccx", "ccy", "ccz", "cswap")
                    length(qubits) == 3 || throw(ArgumentError("$name requires three qubits."))
                    a, b, c = qubits
                    push!(ops, name == "cswap" ? Op("CSWAP", b, a, c) : Op(uppercase(name), a, b, c))
                elseif (!parameterized && name in ("swap", "iswap", "fswap", "syc", "ecr")) ||
                       (parameterized && name in ("rxx", "ryy", "rzz", "rxy", "givens", "fsim", "swapa"))
                    length(qubits) == 2 || throw(ArgumentError("$name requires two qubits."))
                    gate_name = parameterized ? "$name($arguments)" : uppercase(name)
                    push!(ops, Op(gate_name, qubits...))
                else
                    throw(ArgumentError("Unsupported QASM statement: $statement"))
                end
            end
        catch error
            error isa InterruptException && rethrow()
            throw(ArgumentError("$source_name:$line: $(sprint(showerror,error))\n  $statement;"))
        end
    end

    # The vector API carries its register size through the highest used index.
    # Retain a declared idle tail so later compile/run calls allocate every qubit.
    nqubits > get_N(ops) && push!(ops, Op("I", nqubits))
    return ops
end
