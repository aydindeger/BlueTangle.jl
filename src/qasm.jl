###############
# OpenQASM I/O
###############

# ------------ helpers ------------
const _ID = "[A-Za-z_][A-Za-z0-9_]*"
const _QREF = "(" * _ID * ")\\[(\\d+)\\]"                     # qreg[idx]
const _ARGS1 = "\\(\\s*([^(),]+)\\s*\\)"
const _ARGS2 = "\\(\\s*([^(),]+)\\s*,\\s*([^(),]+)\\s*\\)"
const _ARGS3 = "\\(\\s*([^(),]+)\\s*,\\s*([^(),]+)\\s*,\\s*([^(),]+)\\s*\\)"

_clean_upper(name::AbstractString) = uppercase(String(split(name, '(')[1]))
_split_name_args(s::AbstractString) = begin
    if occursin('(', s)
        head, rest = split(s, '('; limit=2)
        return (uppercase(String(head)), String(rstrip(replace(rest, ')' => ""))))
    else
        return (uppercase(String(s)), nothing)
    end
end
_qi0(i::Int) = i - 1  # 1→0
_qi1(s) = parse(Int, s) + 1  # 0→1

"""
    to_qasm(ops::Vector{Op};
            nqubits::Union{Int,Nothing}=nothing,
            qreg::AbstractString="q",
            creg::AbstractString="c") -> String

Emit OpenQASM 2.0 for a vector of `Op`s, including parametrized gates
(RX/RY/RZ, U1/U2/U3, CP, RXX/RYY/RZZ/RXY, GIVENS, FSIM, SWAPA) and controls.
"""
function to_qasm(ops::Vector{Op};
                 nqubits::Union{Int,Nothing}=nothing,
                 qreg::AbstractString="q",
                 creg::AbstractString="c")::String
    isempty(ops) && throw(ArgumentError("ops vector is empty"))

    # --- helpers ---
    _clean_upper(name::AbstractString) = uppercase(String(split(name, '(')[1]))
    _split_name_args(s::AbstractString) = begin
        if occursin('(', s)
            head, rest = split(s, '('; limit=2)
            return (uppercase(String(head)), String(rstrip(replace(rest, ')' => ""))))
        else
            return (uppercase(String(s)), nothing)
        end
    end
    _qi0(i::Int) = i - 1

    # infer N
    function _maxidx(op::Op)
        m = op.qubit
        if op.target_qubit > 0; m = max(m, op.target_qubit); end
        if op.control      > 0; m = max(m, op.control);      end
        m
    end
    N = isnothing(nqubits) ? maximum(_maxidx.(ops)) : nqubits
    N < 1 && throw(ArgumentError("nqubits must be ≥ 1"))

    # collect opaque decls we need to emit
    opaques = Set{String}()

    # op → QASM lines
    function _op_to_qasm(op::Op)::Vector{String}
        q = op.qubit; t = op.target_qubit; c = op.control
        qi = _qi0(q)
        ti = t >= 1 ? _qi0(t) : t
        ci = c >= 1 ? _qi0(c) : c

        nameU = _clean_upper(op.name)
        base, args = _split_name_args(op.name)

        # measurements
        if nameU == "MZ" || nameU == "M(Z)"
            return ["measure $qreg[$qi] -> $creg[$qi];"]
        elseif nameU == "MX" || nameU == "M(X)"
            return ["h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];"]
        elseif nameU == "MY" || nameU == "M(Y)"
            return ["sdg $qreg[$qi];", "h $qreg[$qi];", "measure $qreg[$qi] -> $creg[$qi];"]
        end

        # 1-qubit (incl. parametrized) + optional control
        if t < 0
            if c >= 1
                # ---- Controlled 1q gates ----
                if base == "X";  return ["cx $qreg[$ci],$qreg[$qi];"] end
                if base == "Y";  return ["cy $qreg[$ci],$qreg[$qi];"] end
                if base == "Z";  return ["cz $qreg[$ci],$qreg[$qi];"] end

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
                if     nameU == "I";      return String[] end
                if     nameU == "X";      return ["x $qreg[$qi];"] end
                if     nameU == "Y";      return ["y $qreg[$qi];"] end
                if     nameU == "Z";      return ["z $qreg[$qi];"] end
                if     nameU == "H";      return ["h $qreg[$qi];"] end
                if     nameU == "S";      return ["s $qreg[$qi];"] end
                if     nameU == "SD";     return ["sdg $qreg[$qi];"] end
                if     nameU == "T";      return ["t $qreg[$qi];"] end
                if     nameU == "TD";     return ["tdg $qreg[$qi];"] end
                if     nameU == "SX" || nameU == "XSQRT"; return ["sx $qreg[$qi];"] end
                throw(ArgumentError("Unsupported 1-qubit op: $(op.name)"))
            end
        end

        # 2q / 3q (param and non-param) — unchanged from your last version
        if base == "CX" || base == "CNOT"
            if c >= 1; return ["ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];"] else return ["cx $qreg[$qi],$qreg[$ti];"] end
        elseif base == "CY"
            if c >= 1; return ["ccy $qreg[$qi],$qreg[$ci],$qreg[$ti];"] else return ["cy $qreg[$qi],$qreg[$ti];"] end
        elseif base == "CZ"
            if c >= 1; return ["ccz $qreg[$qi],$qreg[$ci],$qreg[$ti];"] else return ["cz $qreg[$qi],$qreg[$ti];"] end
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
            push!(opaques, "opaque swapa(a) a,b;"); return ["swapa($(args)) $qreg[$qi],$qreg[$ti];"]
        end

        if base == "CCX";   return ["ccx $qreg[$qi],$qreg[$ci],$qreg[$ti];"] end
        if base == "CCY";   return ["ccy $qreg[$qi],$qreg[$ci],$qreg[$ti];"] end
        if base == "CCZ";   return ["ccz $qreg[$qi],$qreg[$ci],$qreg[$ti];"] end
        if base == "CSWAP"; return ["cswap $qreg[$ci],$qreg[$qi],$qreg[$ti];"] end

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


########################
# Robust OpenQASM → Ops
########################

# minimal helpers
const _ID = "[A-Za-z_][A-Za-z0-9_]*"

# 0-based qasm index → 1-based package index, from any token like "q[0]" or "anc[12]"
_q1(tok::AbstractString)::Int = begin
    t = strip(tok)
    # tolerate trailing commas/semicolons/whitespace and any register name
    m = match(r"\[(\d+)\]\s*[;,]?\s*$", t)
    if m === nothing
        throw(ArgumentError("Expected qubit like q[<int>]; got: $tok"))
    end
    parse(Int, m.captures[1]) + 1
end

_qi1(s::AbstractString) = parse(Int, s) + 1   # 0→1

_q2(argstr::AbstractString) = begin
    parts = split(argstr, ',')
    length(parts) == 2 || throw(ArgumentError("Expected two qubits, got: $argstr"))
    (_q1(parts[1]), _q1(parts[2]))
end

function from_qasm(qasm::AbstractString)::Vector{Op}
    ops = Op[]

    # strip // comments and normalize
    src = replace(qasm, r"//.*" => "")
    src = join(split(src, '\n') .|> strip, '\n')
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

        # measure
        if (m = match(r"^measure\s+([^\s]+)\s*->\s*([^\s]+)\s*;", s)) !== nothing
            qi = _q1(m.captures[1]); push!(ops, Op("MZ", qi)); continue
        end

        # unary no-param
        if (m = match(r"^(x|y|z|h|s|sdg|t|tdg|sx)\s+([^\s;]+)\s*;", s)) !== nothing
            kw = m.captures[1]; qi = _q1(m.captures[2])
            oname = kw == "sdg" ? "SD" : kw == "tdg" ? "TD" : kw == "sx" ? "XSQRT" : uppercase(kw)
            push!(ops, Op(oname, qi)); continue
        end

        # unary param: rx/ry/rz, u1/u2/u3
        if (m = match(r"^(rx|ry|rz)\s*\(\s*(.+?)\s*\)\s+([^\s;]+)\s*;", s)) !== nothing
            gate = m.captures[1]; θ = strip(m.captures[2]); qi = _q1(m.captures[3])
            push!(ops, Op("$(gate)($θ)", qi)); continue
        end
        if (m = match(r"^u1\s*\(\s*(.+?)\s*\)\s+([^\s;]+)\s*;", s)) !== nothing
            λ = strip(m.captures[1]); qi = _q1(m.captures[2]); push!(ops, Op("U1($λ)", qi)); continue
        end
        if (m = match(r"^u2\s*\(\s*(.+?)\s*,\s*(.+?)\s*\)\s+([^\s;]+)\s*;", s)) !== nothing
            ϕ = strip(m.captures[1]); λ = strip(m.captures[2]); qi = _q1(m.captures[3])
            push!(ops, Op("U2($ϕ,$λ)", qi)); continue
        end
        if (m = match(r"^u3\s*\(\s*(.+?)\s*,\s*(.+?)\s*,\s*(.+?)\s*\)\s+([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ϕ = strip(m.captures[2]); λ = strip(m.captures[3]); qi = _q1(m.captures[4])
            push!(ops, Op("U3($θ,$ϕ,$λ)", qi)); continue
        end

        # controlled 1q rotations: cu1/cu3 and crx/cry/crz
        if (m = match(r"^(?:cu1|cp)\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            λ = strip(m.captures[1]); ctrl = _q1(m.captures[2]); tgt = _q1(m.captures[3])
            push!(ops, Op("U1($λ)", tgt; control=ctrl)); continue
        end
        if (m = match(r"^cu3\s*\(\s*(.+?)\s*,\s*(.+?)\s*,\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ϕ = strip(m.captures[2]); λ = strip(m.captures[3])
            ctrl = _q1(m.captures[4]); tgt = _q1(m.captures[5])
            if ϕ == "0" && λ == "0"
                push!(ops, Op("ry($θ)", tgt; control=ctrl))
            else
                push!(ops, Op("U3($θ,$ϕ,$λ)", tgt; control=ctrl))
            end
            continue
        end
        if (m = match(r"^crx\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ctrl = _q1(m.captures[2]); tgt = _q1(m.captures[3])
            push!(ops, Op("rx($θ)", tgt; control=ctrl)); continue
        end
        if (m = match(r"^cry\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ctrl = _q1(m.captures[2]); tgt = _q1(m.captures[3])
            push!(ops, Op("ry($θ)", tgt; control=ctrl)); continue
        end
        if (m = match(r"^crz\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ctrl = _q1(m.captures[2]); tgt = _q1(m.captures[3])
            push!(ops, Op("rz($θ)", tgt; control=ctrl)); continue
        end

        # binary non-param (+ map C* to controlled-1q)
        if (m = match(r"^(cx|cnot|cy|cz|swap|iswap|fswap|syc|ecr)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            head = m.captures[1]; a = _q1(m.captures[2]); b = _q1(m.captures[3])
            if head == "cx" || head == "cnot"
                push!(ops, Op("X", b; control=a))
            elseif head == "cy"
                push!(ops, Op("Y", b; control=a))
            elseif head == "cz"
                push!(ops, Op("Z", b; control=a))
            else
                push!(ops, Op(uppercase(head), a, b))
            end
            continue
        end

        # ternary non-param
        if (m = match(r"^(ccx|ccy|ccz|cswap)\s+([^\s,;]+)\s*,\s*([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            head = m.captures[1]
            a = _q1(m.captures[2]); b = _q1(m.captures[3]); c = _q1(m.captures[4])
            if head == "ccx";   push!(ops, Op("CCX", a, b, c))
            elseif head == "ccy"; push!(ops, Op("CCY", a, b, c))
            elseif head == "ccz"; push!(ops, Op("CCZ", a, b, c))
            else                  push!(ops, Op("CSWAP", b, a, c))
            end
            continue
        end

        # binary param: rxx/ryy/rzz/rxy, fsim, givens, swapa
        if (m = match(r"^(rxx|ryy|rzz|rxy)\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            kind = uppercase(m.captures[1]); θ = strip(m.captures[2])
            a = _q1(m.captures[3]); b = _q1(m.captures[4])
            push!(ops, Op("$(kind)($θ)", a, b)); continue
        end
        if (m = match(r"^givens\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); a = _q1(m.captures[2]); b = _q1(m.captures[3])
            push!(ops, Op("GIVENS($θ)", a, b)); continue
        end
        if (m = match(r"^fsim\s*\(\s*(.+?)\s*,\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            θ = strip(m.captures[1]); ϕ = strip(m.captures[2]); a = _q1(m.captures[3]); b = _q1(m.captures[4])
            push!(ops, Op("FSIM($θ,$ϕ)", a, b)); continue
        end
        if (m = match(r"^swapa\s*\(\s*(.+?)\s*\)\s+([^\s,;]+)\s*,\s*([^\s;]+)\s*;", s)) !== nothing
            α = strip(m.captures[1]); a = _q1(m.captures[2]); b = _q1(m.captures[3])
            push!(ops, Op("SWAPA($α)", a, b)); continue
        end

        # if we got here, it’s something we don’t recognize
        throw(ArgumentError("Unsupported or unrecognized statement: $stmt"))
    end

    return ops
end