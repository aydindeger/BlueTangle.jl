###############
# OpenQASM I/O
###############

"""
    to_openqasm(ops::Vector{Op}; nqubits::Union{Int,Nothing}=nothing, qreg::AbstractString="q", creg::AbstractString="c") -> String

Converts a vector of quantum operations (`ops`) into an OpenQASM string representation.

# Arguments
- `ops::Vector{Op}`: A vector containing quantum operations to be converted.
- `nqubits::Union{Int,Nothing}`: (Optional) The number of qubits in the circuit. If `nothing`, the number is inferred from `ops`.
- `qreg::AbstractString`: (Optional) The name of the quantum register. Defaults to `"q"`.
- `creg::AbstractString`: (Optional) The name of the classical register. Defaults to `"c"`.

# Returns
- `String`: The OpenQASM code as a string.

# Throws
- `ArgumentError`: If the `ops` vector is empty.
"""
function to_openqasm(ops::Vector{Op};
                         nqubits::Union{Int,Nothing}=nothing,
                         qreg::AbstractString="q",
                         creg::AbstractString="c")::String
    isempty(ops) && throw(ArgumentError("ops vector is empty"))

    # infer N from all positive indices used
    function _maxidx(op::Op)
        m = op.qubit
        if op.target_qubit > 0; m = max(m, op.target_qubit); end
        if op.control      > 0; m = max(m, op.control);      end
        return m
    end
    N = isnothing(nqubits) ? maximum(_maxidx.(ops)) : nqubits
    N < 1 && throw(ArgumentError("nqubits must be ≥ 1"))

    used_ecr   = any(uppercase(op.name) == "ECR"   for op in ops)
    used_fswap = any(uppercase(op.name) == "FSWAP" for op in ops)

    # Op -> one or more QASM lines
    function _op_to_qasm(op::Op)::Vector{String}
        q = op.qubit; t = op.target_qubit; c = op.control
        qi = q - 1
        ti = t >= 1 ? t - 1 : t
        ci = c >= 1 ? c - 1 : c
        name = uppercase(op.name)

        # Measurements
        if name == "MZ" || name == "M(Z)"
            return ["measure $qreg[$qi] -> $creg[$qi];"]
        elseif name == "MX" || name == "M(X)"
            return ["h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];"]
        elseif name == "MY" || name == "M(Y)"
            return ["sdg $qreg[$qi];", "h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];"]
        end

        # 1-qubit gates (honour optional control)
        if t < 0
            if c >= 1
                if name == "X"
                    return ["cx $qreg[$ci],$qreg[$qi];"]  # control, target
                elseif name == "Y"
                    return ["cy $qreg[$ci],$qreg[$qi];"]
                elseif name == "Z"
                    return ["cz $qreg[$ci],$qreg[$qi];"]
                else
                    throw(ArgumentError("Controlled $(op.name) not supported for export (only X/Y/Z)."))
                end
            else
                if name == "I";                     return String[]               end
                if name == "X";                     return ["x $qreg[$qi];"]      end
                if name == "Y";                     return ["y $qreg[$qi];"]      end
                if name == "Z";                     return ["z $qreg[$qi];"]      end
                if name == "H";                     return ["h $qreg[$qi];"]      end
                if name == "S";                     return ["s $qreg[$qi];"]      end
                if name == "SD";                    return ["sdg $qreg[$qi];"]    end
                if name == "T";                     return ["t $qreg[$qi];"]      end
                if name == "TD";                    return ["tdg $qreg[$qi];"]    end
                if name == "SX" || name == "XSQRT"; return ["sx $qreg[$qi];"]     end
                throw(ArgumentError("Unsupported 1-qubit op: $(op.name)"))
            end
        end

        # 2-/3-qubit and controlled variants by name
        if name == "CX" || name == "CNOT"
            if c >= 1
                return ["ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
            else
                return ["cx $qreg[$qi],$qreg[$ti];"]
            end
        elseif name == "CY"
            if c >= 1
                return ["ccy $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
            else
                return ["cy $qreg[$qi],$qreg[$ti];"]
            end
        elseif name == "CZ"
            if c >= 1
                return ["ccz $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
            else
                return ["cz $qreg[$qi],$qreg[$ti];"]
            end
        elseif name == "SWAP"
            if c >= 1
                return ["cswap $qreg[$ci],$qreg[$qi],$qreg[$ti];"]
            else
                return ["swap $qreg[$qi],$qreg[$ti];"]
            end
        elseif name == "FSWAP"
            return ["fswap $qreg[$qi],$qreg[$ti];"]
        elseif name == "ECR"
            return ["ecr $qreg[$qi],$qreg[$ti];"]
        elseif name == "CCX"
            return ["ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
        elseif name == "CCY"
            return ["ccy $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
        elseif name == "CCZ"
            return ["ccz $qreg[$qi],$qreg[$ci],$qreg[$ti];"]
        elseif name == "CSWAP"
            return ["cswap $qreg[$ci],$qreg[$qi],$qreg[$ti];"]
        else
            throw(ArgumentError("Unsupported multi-qubit op: $(op.name)"))
        end
    end

    lines = String[]
    push!(lines, "OPENQASM 2.0;")
    push!(lines, "include \"qelib1.inc\";")
    if used_ecr;   push!(lines, "opaque ecr a,b;");   end
    if used_fswap; push!(lines, "opaque fswap a,b;"); end
    push!(lines, "qreg $qreg[$N];")
    push!(lines, "creg $creg[$N];")

    for op in ops
        append!(lines, _op_to_qasm(op))
    end
    return join(lines, '\n') * '\n'
end


"""
    from_openqasm(qasm::AbstractString) -> Vector{Op}

Parses an OpenQASM string and returns a vector of quantum operations (`Op`).

# Arguments
- `qasm::AbstractString`: The OpenQASM code as a string.

# Returns
- `Vector{Op}`: A vector containing the parsed quantum operations.
"""
function from_openqasm(qasm::AbstractString)::Vector{Op}
    ops = Op[]

    # strip // comments and normalize whitespace
    src = replace(qasm, r"//.*" => "")
    src = join(split(src, '\n') .|> strip, '\n')

    # helpers (QASM uses 0-based)
    _q1(arg) = begin
        m = match(r"\w+\[(\d+)\]", arg)
        isnothing(m) && throw(ArgumentError("Bad qubit reference: $arg"))
        parse(Int, m.captures[1]) + 1
    end
    _q2(args) = begin
        parts = split(args, ',')
        length(parts) == 2 || throw(ArgumentError("Expected two qubits, got: $args"))
        (_q1(strip(parts[1])), _q1(strip(parts[2])))
    end
    _q3(args) = begin
        parts = split(args, ',')
        length(parts) == 3 || throw(ArgumentError("Expected three qubits, got: $args"))
        (_q1(strip(parts[1])), _q1(strip(parts[2])), _q1(strip(parts[3])))
    end

    stmts = [strip(s)*";" for s in split(src, ';') if !isempty(strip(s))]

    for stmt in stmts
        s = lowercase(stmt)

        # ignore boilerplate
        if startswith(s, "openqasm") || startswith(s, "include") ||
           startswith(s, "qreg")     || startswith(s, "creg")    ||
           startswith(s, "opaque")   || startswith(s, "gate")    ||
           startswith(s, "barrier")
            continue
        end

        # measurements
        if startswith(s, "measure ")
            m = match(r"^measure\s+([^\s]+)\s*->\s*([^\s]+)\s*;", s)
            isnothing(m) && throw(ArgumentError("Bad measure statement: $stmt"))
            qi = _q1(m.captures[1])
            push!(ops, Op("MZ", qi))
            continue
        end

        handled = false

        # unary
        for (qkw, oname) in (
            "x"   => "X",
            "y"   => "Y",
            "z"   => "Z",
            "h"   => "H",
            "s"   => "S",
            "sdg" => "SD",
            "t"   => "T",
            "tdg" => "TD",
            "sx"  => "XSQRT",
        )
            if startswith(s, qkw * " ")
                m = match(r"^\w+\s+([^\s;]+)\s*;", s)
                isnothing(m) && throw(ArgumentError("Bad unary statement: $stmt"))
                qi = _q1(m.captures[1])
                push!(ops, Op(oname, qi))
                handled = true
                break
            end
        end
        handled && continue

        # binary (controlled single-qubit & swaps/ecr/fswap)
        if startswith(s, "cx ") || startswith(s, "cnot ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (ctrl, tgt) = _q2(m.captures[1])
            push!(ops, Op("X", tgt; control=ctrl))
            continue
        elseif startswith(s, "cy ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (ctrl, tgt) = _q2(m.captures[1])
            push!(ops, Op("Y", tgt; control=ctrl))
            continue
        elseif startswith(s, "cz ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (ctrl, tgt) = _q2(m.captures[1])
            push!(ops, Op("Z", tgt; control=ctrl))
            continue
        elseif startswith(s, "swap ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a, b) = _q2(m.captures[1])
            push!(ops, Op("SWAP", a, b))
            continue
        elseif startswith(s, "fswap ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a, b) = _q2(m.captures[1])
            push!(ops, Op("FSWAP", a, b))
            continue
        elseif startswith(s, "ecr ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a, b) = _q2(m.captures[1])
            push!(ops, Op("ECR", a, b))
            continue
        end

        # ternary gates
        if startswith(s, "ccx ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a,b,c) = _q3(m.captures[1])  # controls a,b; target c
            push!(ops, Op("CCX", a, b, c))
            continue
        elseif startswith(s, "ccy ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a,b,c) = _q3(m.captures[1])
            push!(ops, Op("CCY", a, b, c))
            continue
        elseif startswith(s, "ccz ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (a,b,c) = _q3(m.captures[1])
            push!(ops, Op("CCZ", a, b, c))
            continue
        elseif startswith(s, "cswap ")
            m = match(r"^\w+\s+([^\s;]+)\s*;", s)
            (ctrl,q1,q2) = _q3(m.captures[1])  # ctrl, a, b
            push!(ops, Op("CSWAP", q1, ctrl, q2))
            continue
        end

        throw(ArgumentError("Unsupported or unrecognized statement: $stmt"))
    end

    return ops
end