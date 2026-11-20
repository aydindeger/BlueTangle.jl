# Clifford-Augmented Matrix Product State (CAMPS) backend.
#
# Represents an N-qubit pure state as |Ψ⟩ = C |φ⟩, where C is an N-qubit
# Clifford unitary tracked exactly (up to global phase) in O(N²) space by a
# `CliffordTableau`, and |φ⟩ is a disentangled non-Clifford MPS core.

# =============================================================================
# Pauli Strings & Clifford Tableau
# =============================================================================

Base.copy(p::_PauliString) = _PauliString(p.phase, copy(p.x), copy(p.z))
Base.:(==)(a::_PauliString, b::_PauliString) = a.phase == b.phase && a.x == b.x && a.z == b.z

function _pauli_single(n::Int, q::Int, kind::Symbol; phase::Int=1)
    x = falses(n)
    z = falses(n)
    if kind === :X
        x[q] = true
    elseif kind === :Z
        z[q] = true
    elseif kind === :Y
        x[q] = true
        z[q] = true
    elseif kind !== :I
        throw(ArgumentError("unknown Pauli kind: $kind"))
    end
    return _PauliString(Int8(phase), x, z)
end

function _pauli_weight(p::_PauliString)
    w = 0
    @inbounds for q in eachindex(p.x)
        (p.x[q] || p.z[q]) && (w += 1)
    end
    return w
end

function _pauli_to_string(p::_PauliString)
    chars = Vector{Char}(undef, length(p.x))
    @inbounds for q in eachindex(p.x)
        xq, zq = p.x[q], p.z[q]
        chars[q] = (xq && zq) ? 'Y' : (xq ? 'X' : (zq ? 'Z' : 'I'))
    end
    prefix = p.phase == 1 ? "+" : "-"
    return prefix * String(chars)
end

# Multiply two Hermitian Pauli strings: A * B = (i^ipow) * C with ipow ∈ {0, 1}.
function _pauli_mul_raw(a::_PauliString, b::_PauliString)
    n = length(a.x)
    x = a.x .⊻ b.x
    z = a.z .⊻ b.z
    ipow = 0
    @inbounds for q in 1:n
        ax, az = a.x[q], a.z[q]
        bx, bz = b.x[q], b.z[q]
        if ax && !az          # A = X
            if !bx && bz      # X * Z = -i Y
                ipow += 3
            elseif bx && bz   # X * Y = +i Z
                ipow += 1
            end
        elseif !ax && az      # A = Z
            if bx && !bz      # Z * X = +i Y
                ipow += 1
            elseif bx && bz   # Z * Y = -i X
                ipow += 3
            end
        elseif ax && az       # A = Y
            if bx && !bz      # Y * X = -i Z
                ipow += 3
            elseif !bx && bz  # Y * Z = +i X
                ipow += 1
            end
        end
    end
    ipow_mod = mod(ipow, 4)
    ph = Int(a.phase) * Int(b.phase)
    if ipow_mod == 2
        ph = -ph
        ipow_mod = 0
    elseif ipow_mod == 3
        ph = -ph
        ipow_mod = 1
    end
    return ipow_mod, _PauliString(Int8(ph), x, z)
end

function _pauli_mul_commuting(a::_PauliString, b::_PauliString)
    ipow, c = _pauli_mul_raw(a, b)
    ipow == 0 || error("expected commuting Pauli strings")
    return c
end

# Compute (imag_sign * i) * A * B for anticommuting Hermitian Pauli strings.
function _pauli_mul_anticommuting(a::_PauliString, b::_PauliString; imag_sign::Int=1)
    ipow, c = _pauli_mul_raw(a, b)
    ipow == 1 || error("expected anticommuting Pauli strings")
    # A * B = i * c => (imag_sign * i) * (i * c) = (-imag_sign) * c
    return _PauliString(Int8(-imag_sign * Int(c.phase)), c.x, c.z)
end

function CliffordTableau(n::Int)
    n > 0 || throw(ArgumentError("number of qubits must be positive"))
    return CliffordTableau(
        n,
        [_pauli_single(n, j, :X) for j in 1:n],
        [_pauli_single(n, j, :Z) for j in 1:n],
    )
end

Base.copy(tab::CliffordTableau) = CliffordTableau(tab.n, copy.(tab.px), copy.(tab.pz))
Base.:(==)(a::CliffordTableau, b::CliffordTableau) = a.n == b.n && a.px == b.px && a.pz == b.pz

function Base.show(io::IO, tab::CliffordTableau)
    print(io, "CliffordTableau(N=$(tab.n))")
end

function Base.show(io::IO, ::MIME"text/plain", tab::CliffordTableau)
    println(io, "CliffordTableau(N=$(tab.n))")
    limit = min(tab.n, 8)
    for j in 1:limit
        println(io, "  C† X_$j C = ", _pauli_to_string(tab.px[j]),
                "    C† Z_$j C = ", _pauli_to_string(tab.pz[j]))
    end
    tab.n > limit && println(io, "  ⋮ ($(tab.n - limit) more qubits)")
end

function CAMPSState(n::Int)
    n > 0 || throw(ArgumentError("number of qubits must be positive"))
    return CAMPSState(CliffordTableau(n), zero_state(N_MPS(n)))
end

function CAMPSState(psi::it.MPS)
    n = length(psi)
    n > 0 || throw(ArgumentError("MPS must have at least one site"))
    return CAMPSState(CliffordTableau(n), copy(psi))
end

Base.copy(st::CAMPSState) = CAMPSState(copy(st.tableau), copy(st.mps))

function Base.show(io::IO, st::CAMPSState)
    print(io, "CAMPSState(N=$(get_N(st)), maxlinkdim=$(maxlinkdim(st)))")
end

function _tableau_is_identity(tab::CliffordTableau)
    @inbounds for j in 1:tab.n
        px, pz = tab.px[j], tab.pz[j]
        (px.phase == 1 && pz.phase == 1) || return false
        for q in 1:tab.n
            expected = (q == j)
            (px.x[q] == expected && !px.z[q]) || return false
            (!pz.x[q] && pz.z[q] == expected) || return false
        end
    end
    return true
end

# Left-multiply C ← U * C by an elementary Clifford gate U.
# Updates rows C† P_j C ← C† (U† P_j U) C.
function _tableau_left_gate!(tab::CliffordTableau, name::AbstractString, q::Int, t::Int=-1)
    u = uppercase(String(name))
    if u == "I" || u == "CI"
        return tab
    elseif u == "H"
        tab.px[q], tab.pz[q] = tab.pz[q], tab.px[q]
    elseif u == "S"
        # S† X S = -Y = -i X Z
        tab.px[q] = _pauli_mul_anticommuting(tab.px[q], tab.pz[q]; imag_sign=-1)
    elseif u == "SD" || u == "SDG"
        # S X S† = +Y = +i X Z
        tab.px[q] = _pauli_mul_anticommuting(tab.px[q], tab.pz[q]; imag_sign=+1)
    elseif u == "SX" || u == "XSQRT"
        # SX = H * S * H => SX† Z SX = +Y = +i X Z, SX† X SX = X
        tab.pz[q] = _pauli_mul_anticommuting(tab.px[q], tab.pz[q]; imag_sign=+1)
    elseif u in ("SXDG", "SXD", "XSQRTD", "XSQRTDG")
        # SXDG† Z SXDG = -Y = -i X Z, SXDG† X SXDG = X
        tab.pz[q] = _pauli_mul_anticommuting(tab.px[q], tab.pz[q]; imag_sign=-1)
    elseif u == "HSP"
        # HSP = H * S†
        _tableau_left_gate!(tab, "SD", q)
        _tableau_left_gate!(tab, "H", q)
    elseif u == "HY"
        # HY = H * S
        _tableau_left_gate!(tab, "S", q)
        _tableau_left_gate!(tab, "H", q)
    elseif u == "H2"
        _tableau_left_gate!(tab, "H", q)
        _tableau_left_gate!(tab, "H", t)
    elseif u == "X"
        tab.pz[q].phase = -tab.pz[q].phase
    elseif u == "Z"
        tab.px[q].phase = -tab.px[q].phase
    elseif u == "Y"
        tab.px[q].phase = -tab.px[q].phase
        tab.pz[q].phase = -tab.pz[q].phase
    elseif u == "CX" || u == "CNOT"
        tab.px[q] = _pauli_mul_commuting(tab.px[q], tab.px[t])
        tab.pz[t] = _pauli_mul_commuting(tab.pz[q], tab.pz[t])
    elseif u == "CZ"
        tab.px[q] = _pauli_mul_commuting(tab.px[q], tab.pz[t])
        tab.px[t] = _pauli_mul_commuting(tab.pz[q], tab.px[t])
    elseif u == "CY"
        # CY(q, t) = S_t * CX(q, t) * SD_t
        _tableau_left_gate!(tab, "SD", t)
        _tableau_left_gate!(tab, "CX", q, t)
        _tableau_left_gate!(tab, "S", t)
    elseif u == "SWAP"
        tab.px[q], tab.px[t] = tab.px[t], tab.px[q]
        tab.pz[q], tab.pz[t] = tab.pz[t], tab.pz[q]
    elseif u == "ISWAP"
        # ISWAP = SWAP * CZ * S_q * S_t
        _tableau_left_gate!(tab, "S", q)
        _tableau_left_gate!(tab, "S", t)
        _tableau_left_gate!(tab, "CZ", q, t)
        _tableau_left_gate!(tab, "SWAP", q, t)
    elseif u == "FSWAP"
        # FSWAP = SWAP * CZ
        _tableau_left_gate!(tab, "CZ", q, t)
        _tableau_left_gate!(tab, "SWAP", q, t)
    elseif u == "ECR"
        # ECR = (IX - XY)/√2 (Hermitian Clifford):
        # ECR† X_q ECR = X_q, ECR† Z_q ECR = Y_q Z_t,
        # ECR† X_t ECR = -X_q Y_t, ECR† Z_t ECR = -Z_t.
        pxq, pzq, pxt, pzt = tab.px[q], tab.pz[q], tab.px[t], tab.pz[t]
        yq = _pauli_mul_anticommuting(pxq, pzq; imag_sign=+1)
        yt = _pauli_mul_anticommuting(pxt, pzt; imag_sign=+1)
        new_pzq = _pauli_mul_commuting(yq, pzt)
        new_pxt = _pauli_mul_commuting(pxq, yt)
        new_pxt.phase = -new_pxt.phase
        new_pzt = copy(pzt)
        new_pzt.phase = -new_pzt.phase
        tab.pz[q] = new_pzq
        tab.px[t] = new_pxt
        tab.pz[t] = new_pzt
    else
        throw(ArgumentError("unsupported Clifford gate in _tableau_left_gate!: $name"))
    end
    return tab
end

# Left-multiply C ← P * C for an arbitrary Pauli string P.
function _tableau_left_pauli!(tab::CliffordTableau, p::_PauliString)
    @inbounds for j in 1:tab.n
        p.z[j] && (tab.px[j].phase = -tab.px[j].phase)
        p.x[j] && (tab.pz[j].phase = -tab.pz[j].phase)
    end
    return tab
end

# Left-multiply C ← exp(-i * sign * (π/4) * P) * C where sign ∈ {+1, -1}.
function _tableau_left_pauli_phase!(tab::CliffordTableau, p::_PauliString, sign::Int)
    p_tilde = _tableau_pullback(tab, p)
    @inbounds for j in 1:tab.n
        if p.z[j]
            tab.px[j] = _pauli_mul_anticommuting(p_tilde, tab.px[j]; imag_sign=sign)
        end
        if p.x[j]
            tab.pz[j] = _pauli_mul_anticommuting(p_tilde, tab.pz[j]; imag_sign=sign)
        end
    end
    return tab
end

# Conjugate a Pauli string in place: P ← D P D† for elementary Clifford gate D.
function _pauli_conj_gate!(p::_PauliString, name::AbstractString, q::Int, t::Int=-1)
    u = uppercase(String(name))
    xq, zq = p.x[q], p.z[q]
    if u == "I"
        return p
    elseif u == "H"
        p.x[q], p.z[q] = zq, xq
        (xq && zq) && (p.phase = -p.phase)
    elseif u == "S"
        (xq && zq) && (p.phase = -p.phase)
        p.z[q] = zq ⊻ xq
    elseif u == "SD" || u == "SDG"
        (xq && !zq) && (p.phase = -p.phase)
        p.z[q] = zq ⊻ xq
    elseif u == "X"
        zq && (p.phase = -p.phase)
    elseif u == "Z"
        xq && (p.phase = -p.phase)
    elseif u == "Y"
        (xq ⊻ zq) && (p.phase = -p.phase)
    elseif u == "CX" || u == "CNOT"
        xt, zt = p.x[t], p.z[t]
        if xq && zt && !(xt ⊻ zq)
            p.phase = -p.phase
        end
        p.x[t] = xt ⊻ xq
        p.z[q] = zq ⊻ zt
    elseif u == "CZ"
        _pauli_conj_gate!(p, "H", t)
        _pauli_conj_gate!(p, "CX", q, t)
        _pauli_conj_gate!(p, "H", t)
    elseif u == "CY"
        _pauli_conj_gate!(p, "SD", t)
        _pauli_conj_gate!(p, "CX", q, t)
        _pauli_conj_gate!(p, "S", t)
    elseif u == "SWAP"
        p.x[q], p.x[t] = p.x[t], p.x[q]
        p.z[q], p.z[t] = p.z[t], p.z[q]
    else
        throw(ArgumentError("unsupported Clifford gate in _pauli_conj_gate!: $name"))
    end
    return p
end

# Right-multiply C ← C * D† (paired with |φ⟩ ← D |φ⟩ so C |φ⟩ is invariant).
function _tableau_right_dag_gate!(tab::CliffordTableau, name::AbstractString, q::Int, t::Int=-1)
    @inbounds for j in 1:tab.n
        _pauli_conj_gate!(tab.px[j], name, q, t)
        _pauli_conj_gate!(tab.pz[j], name, q, t)
    end
    return tab
end

# Compute C† P C for any Hermitian Pauli string P.
function _tableau_pullback(tab::CliffordTableau, p::_PauliString)
    acc = _PauliString(p.phase, falses(tab.n), falses(tab.n))
    @inbounds for q in 1:tab.n
        xq, zq = p.x[q], p.z[q]
        if xq && !zq
            acc = _pauli_mul_commuting(acc, tab.px[q])
        elseif !xq && zq
            acc = _pauli_mul_commuting(acc, tab.pz[q])
        elseif xq && zq
            py = _pauli_mul_anticommuting(tab.px[q], tab.pz[q]; imag_sign=+1)
            acc = _pauli_mul_commuting(acc, py)
        end
    end
    return acc
end

function _pauli_to_ops(p::_PauliString)
    ops = Op[]
    @inbounds for q in eachindex(p.x)
        xq, zq = p.x[q], p.z[q]
        if xq && zq
            push!(ops, Op("Y", q; noisy=false))
        elseif xq
            push!(ops, Op("X", q; noisy=false))
        elseif zq
            push!(ops, Op("Z", q; noisy=false))
        end
    end
    return ops
end

"""
    clifford_ops(tableau::CliffordTableau) -> Vector{Op}
    clifford_ops(state::CAMPSState) -> Vector{Op}

Synthesize a sequence of elementary Clifford `Op`s (`H`, `S`, `SD`, `X`, `Y`,
`Z`, `CX`, `CZ`, `CY`, `SWAP`, all with `noisy=false`) realizing the Clifford
unitary `C` tracked by `tableau` up to a single global phase, using `O(N²)`
symplectic Gaussian elimination.
"""
function clifford_ops(tab::CliffordTableau)
    t = copy(tab)
    n = t.n
    ops = Op[]

    function emit!(name::String, q::Int, target::Int=-1)
        _tableau_right_dag_gate!(t, name, q, target)
        push!(ops, target < 0 ? Op(name, q; noisy=false) : Op(name, q, target; noisy=false))
    end

    for k in 1:n
        # Step 1: reduce t.pz[k] to ±Z_k using only gates on qubits k:n.
        pz = t.pz[k]
        support = Int[]
        for q in k:n
            if pz.x[q] || pz.z[q]
                push!(support, q)
            end
        end
        isempty(support) && error("corrupted CliffordTableau: identity Z row at qubit $k")
        for q in support
            xq, zq = pz.x[q], pz.z[q]
            if xq && zq
                emit!("S", q)
                emit!("H", q)
            elseif xq
                emit!("H", q)
            end
        end
        pivot = first(support)
        for idx in 2:length(support)
            emit!("CX", support[idx], pivot)
        end
        if pivot != k
            emit!("SWAP", k, pivot)
        end

        # Step 2: reduce t.px[k] to ±X_k while keeping t.pz[k] = ±Z_k.
        px = t.px[k]
        if px.x[k] && px.z[k]
            emit!("SD", k)
        end
        for q in (k + 1):n
            xq, zq = px.x[q], px.z[q]
            if xq && zq
                emit!("CY", k, q)
            elseif xq
                emit!("CX", k, q)
            elseif zq
                emit!("CZ", k, q)
            end
        end

        # Step 3: fix signs of t.px[k] and t.pz[k].
        neg_x = (t.px[k].phase == -1)
        neg_z = (t.pz[k].phase == -1)
        if neg_x && neg_z
            emit!("Y", k)
        elseif neg_x
            emit!("Z", k)
        elseif neg_z
            emit!("X", k)
        end
    end

    return ops
end

clifford_ops(st::CAMPSState) = clifford_ops(st.tableau)

# Compute the Clifford tableau representing C_a† * C_b.
function _tableau_relative(tab_a::CliffordTableau, tab_b::CliffordTableau)
    tab_a.n == tab_b.n || throw(DimensionMismatch("tableaux must have the same number of qubits"))
    rel = copy(tab_b)
    ops_a = clifford_ops(tab_a)
    for op in Iterators.reverse(ops_a)
        nm = uppercase(op.name)
        inv_nm = nm == "S" ? "SD" : (nm == "SD" || nm == "SDG") ? "S" : nm
        _tableau_left_gate!(rel, inv_nm, op.qubit, op.target_qubit)
    end
    return rel
end

# =============================================================================
# 2-Qubit Clifford Disentangler Cosets (C₁^{⊗2} \ C₂ : 19 cosets)
# =============================================================================

struct _CAMPSDisentangler2Q
    gates::Vector{Tuple{String,Int,Int}}
    mat::Matrix{ComplexF64}
end

function _camps_2q_elem_mat(name::String, q::Int, t::Int=-1)
    if t < 0
        g = Matrix{ComplexF64}(gates(name))
        I2 = Matrix{ComplexF64}(la.I, 2, 2)
        return q == 1 ? kron(g, I2) : kron(I2, g)
    elseif name == "CX"
        return q == 1 ? Matrix{ComplexF64}(gate.CX) :
               Matrix{ComplexF64}(gate.SWAP * gate.CX * gate.SWAP)
    elseif name == "SWAP"
        return Matrix{ComplexF64}(gate.SWAP)
    else
        error("unexpected disentangler gate: $name")
    end
end

function _build_camps_2q_cosets()
    c_reps = [
        Tuple{String,Int,Int}[],
        [("H", 1, -1)],
        [("S", 1, -1), ("H", 1, -1)],
    ]
    t_reps = [
        Tuple{String,Int,Int}[],
        [("H", 2, -1)],
        [("S", 2, -1)],
    ]

    cosets = _CAMPSDisentangler2Q[]
    function add_seq(seq::Vector{Tuple{String,Int,Int}})
        M = Matrix{ComplexF64}(la.I, 4, 4)
        for (nm, q, t) in seq
            M = _camps_2q_elem_mat(nm, q, t) * M
        end
        push!(cosets, _CAMPSDisentangler2Q(seq, M))
    end

    # Identity (1 coset)
    add_seq(Tuple{String,Int,Int}[])
    # CNOT class (9 cosets)
    for r1 in c_reps, r2 in t_reps
        add_seq(vcat(r1, r2, [("CX", 1, 2)]))
    end
    # Double-CNOT / iSWAP class (9 cosets)
    for r1 in c_reps, r2 in t_reps
        add_seq(vcat(r1, r2, [("CX", 1, 2), ("CX", 2, 1)]))
    end
    return cosets
end

const _CAMPS_COSETS_2Q = _build_camps_2q_cosets()

function _camps_vn_entropy(svals::Vector{Float64})
    p = svals .^ 2
    Z = sum(p)
    Z <= 0 && return 0.0
    S = 0.0
    invZ = 1.0 / Z
    @inbounds for pk in p
        q = pk * invZ
        q > 1e-15 && (S -= q * log(q))
    end
    return S
end

function _camps_trunc_kwargs(maxdim::Union{Nothing,Int}, cutoff::Real; scale::Int=1)
    return maxdim === nothing ? (; cutoff=cutoff) : (; maxdim=scale * maxdim, cutoff=cutoff)
end

function _camps_disentangle_sweep!(st::CAMPSState;
                                   maxdim::Union{Nothing,Int}=nothing,
                                   cutoff::Real=1e-12,
                                   nsweeps::Int=2,
                                   tol_improve::Real=1e-9)
    psi = st.mps
    n = length(psi)
    (n <= 1 || nsweeps <= 0) && return false
    sites = it.siteinds(psi)
    trunc_kw = _camps_trunc_kwargs(maxdim, cutoff)
    improved_any = false

    for sweep in 1:nsweeps
        bonds = isodd(sweep) ? (1:(n - 1)) : ((n - 1):-1:1)
        sweep_improved = false
        for b in bonds
            it.orthogonalize!(psi, b)
            wf = psi[b] * psi[b + 1]
            linds = b == 1 ? (sites[b],) : (it.linkind(psi, b - 1), sites[b])

            U0, S0, _ = it.svd(wf, linds; cutoff=cutoff)
            sv0 = [S0[k, k] for k in 1:it.dim(it.commonind(U0, S0))]
            best_S = _camps_vn_entropy(sv0)
            best_idx = 1
            best_wf = wf

            if length(sv0) > 1 || best_S > 1e-12
                for idx in 2:length(_CAMPS_COSETS_2Q)
                    D = _CAMPS_COSETS_2Q[idx]
                    op_t = it.op(D.mat, sites[b + 1], sites[b])
                    wf_try = it.noprime(op_t * wf)
                    Ut, St, _ = it.svd(wf_try, linds; cutoff=cutoff)
                    svt = [St[k, k] for k in 1:it.dim(it.commonind(Ut, St))]
                    St_ent = _camps_vn_entropy(svt)
                    if St_ent < best_S - tol_improve
                        best_S = St_ent
                        best_idx = idx
                        best_wf = wf_try
                        best_S <= 1e-12 && break
                    end
                end
            end

            if best_idx != 1
                D = _CAMPS_COSETS_2Q[best_idx]
                for (nm, lq, lt) in D.gates
                    gq = lq == 1 ? b : b + 1
                    gt = lt < 0 ? -1 : (lt == 1 ? b : b + 1)
                    _tableau_right_dag_gate!(st.tableau, nm, gq, gt)
                end
                sweep_improved = true
                improved_any = true
            end

            U, S, V = it.svd(best_wf, linds; trunc_kw...)
            if isodd(sweep)
                psi[b] = U
                psi[b + 1] = S * V
            else
                psi[b] = U * S
                psi[b + 1] = V
            end
        end
        sweep_improved || break
    end
    st.mps = la.normalize(psi)
    return improved_any
end

# Identify sites q where |φ⟩ is currently in the product state |0⟩_q.
function _camps_free_zero_sites(psi::it.MPS; tol::Real=1e-11)
    z_exps = real.(it.expect(psi, "Z"))
    zeros_q = Int[]
    @inbounds for q in eachindex(z_exps)
        z_exps[q] > 1.0 - tol && push!(zeros_q, q)
    end
    return zeros_q
end

# Collapse P̃ using untouched |0⟩ sites in |φ⟩ without changing C|φ⟩.
# Returns `pivot` (the free-zero site now carrying ±X_pivot) or `0` if none.
function _camps_reduce_pauli_via_free_zeros!(tab::CliffordTableau, psi::it.MPS, p_tilde::_PauliString)
    free_zeros = _camps_free_zero_sites(psi)
    isempty(free_zeros) && return 0

    pivot = 0
    for q in free_zeros
        if p_tilde.x[q]
            pivot = q
            break
        end
    end

    if pivot != 0
        if p_tilde.z[pivot]
            _tableau_right_dag_gate!(tab, "SD", pivot)
            _pauli_conj_gate!(p_tilde, "SD", pivot)
        end
        for a in 1:tab.n
            a == pivot && continue
            xa, za = p_tilde.x[a], p_tilde.z[a]
            if xa || za
                gname = (xa && za) ? "CY" : (xa ? "CX" : "CZ")
                _tableau_right_dag_gate!(tab, gname, pivot, a)
                _pauli_conj_gate!(p_tilde, gname, pivot, a)
            end
        end
        return pivot
    end

    # P̃ has only I or Z on all free_zeros; drop Z factors since Z_q |0⟩_q = |0⟩_q.
    for q in free_zeros
        if !p_tilde.x[q] && p_tilde.z[q]
            p_tilde.z[q] = false
        end
    end
    return 0
end

# =============================================================================
# Pauli Rotations & Unitary Gate Decompositions
# =============================================================================

# Apply exp(-i * (θ/2) * P) to st = C|φ⟩.
function _camps_apply_pauli_string_rot!(st::CAMPSState, p::_PauliString, theta::Real;
                                        maxdim::Union{Nothing,Int}=nothing,
                                        cutoff::Real=1e-12,
                                        disentangle::Bool=true,
                                        nsweeps::Int=2)
    # Reduce θ modulo 2π (up to global phase) to detect exact Clifford angles.
    eff = Float64(p.phase) * Float64(theta)
    half_pi_units = round(eff / (π / 2))
    if abs(eff - half_pi_units * (π / 2)) < 1e-12
        k = mod(Int(half_pi_units), 4)
        if k == 0
            return st
        elseif k == 2
            _tableau_left_pauli!(st.tableau, p)
            return st
        elseif k == 1
            # exp(-i (π/4) P_unsigned)
            p_unsigned = _PauliString(Int8(1), p.x, p.z)
            _tableau_left_pauli_phase!(st.tableau, p_unsigned, +1)
            return st
        else # k == 3 => -π/2 => exp(+i (π/4) P_unsigned)
            p_unsigned = _PauliString(Int8(1), p.x, p.z)
            _tableau_left_pauli_phase!(st.tableau, p_unsigned, -1)
            return st
        end
    end

    p_tilde = _tableau_pullback(st.tableau, p)
    if disentangle
        _camps_reduce_pauli_via_free_zeros!(st.tableau, st.mps, p_tilde)
    end
    pauli_ops = _pauli_to_ops(p_tilde)

    if isempty(pauli_ops)
        return st
    elseif length(pauli_ops) == 1
        op1 = only(pauli_ops)
        rot_angle = Int(p_tilde.phase) * Float64(theta)
        rot_mat = op1.name == "X" ? _RX(rot_angle) :
                  op1.name == "Y" ? _RY(rot_angle) : _RZ(rot_angle)
        sites = it.siteinds(st.mps)
        st.mps = la.normalize(it.apply(it.op(rot_mat, sites[op1.qubit]), st.mps; cutoff=cutoff))
        return st
    end

    c = cos(theta / 2)
    s = -im * sin(theta / 2) * Int(p_tilde.phase)
    p_phi = _camps_apply_pauli_ops_mps(st.mps, pauli_ops)
    combined = +(c * st.mps, s * p_phi; alg="directsum")

    if disentangle
        pre_kw = _camps_trunc_kwargs(maxdim, cutoff; scale=2)
        it.truncate!(combined; pre_kw...)
        st.mps = la.normalize(combined)
        _camps_disentangle_sweep!(st; maxdim=maxdim, cutoff=cutoff, nsweeps=nsweeps)
    else
        trunc_kw = _camps_trunc_kwargs(maxdim, cutoff)
        it.truncate!(combined; trunc_kw...)
        st.mps = la.normalize(combined)
    end
    return st
end

function _camps_apply_pauli_ops_mps(psi::it.MPS, pauli_ops::Vector{Op})
    isempty(pauli_ops) && return copy(psi)
    sites = it.siteinds(psi)
    out = copy(psi)
    for op in pauli_ops
        out = it.apply(it.op(op.mat, sites[op.qubit]), out; cutoff=0.0)
    end
    return out
end

function _camps_apply_axis_rot!(st::CAMPSState, axis::Symbol, q::Int, theta::Real; kwargs...)
    return _camps_apply_pauli_string_rot!(st, _pauli_single(st.tableau.n, q, axis), theta; kwargs...)
end

# Apply an arbitrary 2×2 unitary matrix U on qubit q.
# Detects single-axis Z, X, and Y rotations first (1 rotation instead of 3),
# and routes all 24 single-qubit Cliffords through the zero-MPS-cost tableau path.
function _camps_apply_1q_unitary!(st::CAMPSState, U::AbstractMatrix, q::Int; kwargs...)
    M = Matrix{ComplexF64}(U)
    isunitary(M) || throw(ArgumentError("CAMPS single-qubit gate matrix must be unitary"))

    # Diagonal => pure Z rotation
    if abs(M[1, 2]) < 1e-12 && abs(M[2, 1]) < 1e-12
        theta = angle(M[2, 2] / M[1, 1])
        return _camps_apply_axis_rot!(st, :Z, q, theta; kwargs...)
    end

    # Diagonal in X basis (H M H diagonal) => pure X rotation
    Mx = gate.H * M * gate.H
    if abs(Mx[1, 2]) < 1e-12 && abs(Mx[2, 1]) < 1e-12
        theta = angle(Mx[2, 2] / Mx[1, 1])
        return _camps_apply_axis_rot!(st, :X, q, theta; kwargs...)
    end

    # Diagonal in Y basis (HSP M HSP' diagonal, since HSP Y HSP' = Z) => pure Y rotation
    My = gate.HSP * M * gate.HSP'
    if abs(My[1, 2]) < 1e-12 && abs(My[2, 1]) < 1e-12
        theta = angle(My[2, 2] / My[1, 1])
        return _camps_apply_axis_rot!(st, :Y, q, theta; kwargs...)
    end

    # General 1-qubit unitary: U = exp(iα) * RZ(β) * RY(γ) * RZ(δ)
    _, β, γ, δ = zyz_decomposition(M)
    abs(δ) > 1e-12 && _camps_apply_axis_rot!(st, :Z, q, δ; kwargs...)
    abs(γ) > 1e-12 && _camps_apply_axis_rot!(st, :Y, q, γ; kwargs...)
    abs(β) > 1e-12 && _camps_apply_axis_rot!(st, :Z, q, β; kwargs...)
    return st
end

function _pauli_two(n::Int, q::Int, k1::Symbol, t::Int, k2::Symbol)
    return _pauli_mul_commuting(_pauli_single(n, q, k1), _pauli_single(n, t, k2))
end

# Exact Kraus-Cirac magic-basis KAK decomposition for any U ∈ U(4):
#   U = e^{iα} (A₁ ⊗ A₂) exp(-i(θx/2 XX + θy/2 YY + θz/2 ZZ)) (B₁ ⊗ B₂).
const _KAK_MAGIC_Q = (1.0 / sqrt(2.0)) * ComplexF64[
    1  0   0  im;
    0  im  1  0;
    0  im -1  0;
    1  0   0 -im
]

function _kak_simultaneous_real_eigen(R::Matrix{Float64}, Imat::Matrix{Float64})
    F = la.eigen(la.Symmetric(R))
    vals = F.values
    V = copy(F.vectors)
    # Break any degenerate eigenspace of R using the commuting real symmetric matrix Imat.
    i = 1
    while i <= 4
        j = i
        while j < 4 && abs(vals[j + 1] - vals[i]) < 1e-8
            j += 1
        end
        if j > i
            subV = V[:, i:j]
            subI = la.Symmetric(subV' * Imat * subV)
            V[:, i:j] = subV * la.eigen(subI).vectors
        end
        i = j + 1
    end
    if la.det(V) < 0
        V[:, 1] .= -V[:, 1]
    end
    return V
end

function _kak_decompose_2q(U::AbstractMatrix)
    M0 = Matrix{ComplexF64}(U)
    isunitary(M0) || throw(ArgumentError("CAMPS two-qubit gate matrix must be unitary"))
    # Normalize to SU(4)
    det_phase = angle(la.det(M0)) / 4
    Usu = M0 * cis(-det_phase)

    M = _KAK_MAGIC_Q' * Usu * _KAK_MAGIC_Q
    W = transpose(M) * M
    O = _kak_simultaneous_real_eigen(real.(W), imag.(W))
    Ddiag = la.diag(transpose(O) * W * O)
    λ = angle.(Ddiag) ./ 2
    Dhalf = la.Diagonal(cis.(λ))

    K1c = M * O * conj(Dhalf)
    K1 = real.(K1c)
    if la.det(K1) < 0
        K1[:, 1] .= -K1[:, 1]
        λ[1] += π
    end

    Amat = _KAK_MAGIC_Q * K1 * _KAK_MAGIC_Q'
    Bmat = _KAK_MAGIC_Q * transpose(O) * _KAK_MAGIC_Q'
    A1, A2 = kronecker_decomposition(Amat)
    B1, B2 = kronecker_decomposition(Bmat)

    θx = (λ[3] + λ[4] - λ[1] - λ[2]) / 2
    θy = (λ[1] + λ[3] - λ[2] - λ[4]) / 2
    θz = (λ[2] + λ[3] - λ[1] - λ[4]) / 2
    return A1, A2, θx, θy, θz, B1, B2
end

function _camps_apply_cp!(st::CAMPSState, q::Int, t::Int, λ::Real; kwargs...)
    n = st.tableau.n
    _camps_apply_axis_rot!(st, :Z, q, λ / 2; kwargs...)
    _camps_apply_axis_rot!(st, :Z, t, λ / 2; kwargs...)
    _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Z, t, :Z), -λ / 2; kwargs...)
    return st
end

function _camps_apply_2q_unitary!(st::CAMPSState, mat::AbstractMatrix, q::Int, t::Int;
                                  name::AbstractString="", kwargs...)
    n = st.tableau.n
    up = uppercase(_clean_name(String(name)))

    # Fast paths for built-in 2-qubit gates when matrix matches the built-in definition.
    if up in ("CI", "CX", "CNOT", "CZ", "CY", "SWAP", "ISWAP", "FSWAP", "ECR", "H2") &&
       isapprox(mat, gates(up); atol=1e-12, rtol=1e-12)
        return _tableau_left_gate!(st.tableau, up, q, t)
    elseif up == "CS" && isapprox(mat, gate.CS; atol=1e-12, rtol=1e-12)
        return _camps_apply_cp!(st, q, t, π / 2; kwargs...)
    elseif up == "CSD" && isapprox(mat, gate.CSD; atol=1e-12, rtol=1e-12)
        return _camps_apply_cp!(st, q, t, -π / 2; kwargs...)
    elseif up == "CT" && isapprox(mat, gate.CT; atol=1e-12, rtol=1e-12)
        return _camps_apply_cp!(st, q, t, π / 4; kwargs...)
    elseif up == "CTD" && isapprox(mat, gate.CTD; atol=1e-12, rtol=1e-12)
        return _camps_apply_cp!(st, q, t, -π / 4; kwargs...)
    elseif up == "CH" && isapprox(mat, gate.CH; atol=1e-12, rtol=1e-12)
        _camps_apply_axis_rot!(st, :Y, t, -π / 4; kwargs...)
        _tableau_left_gate!(st.tableau, "CZ", q, t)
        _camps_apply_axis_rot!(st, :Y, t, π / 4; kwargs...)
        return st
    elseif up == "SYC" && isapprox(mat, gate.SYC; atol=1e-12, rtol=1e-12)
        _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), π / 2; kwargs...)
        _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), π / 2; kwargs...)
        return _camps_apply_cp!(st, q, t, -π / 6; kwargs...)
    elseif up in ("CP", "RXX", "RYY", "RZZ", "RXY", "GIVENS", "FSIM", "SWAPA") &&
           _name_with_arg_bool(String(name)) && isapprox(mat, gates(String(name)); atol=1e-12, rtol=1e-12)
        args = _get_name_and_arguments(String(name))[3]
        if up == "CP"
            return _camps_apply_cp!(st, q, t, args[1]; kwargs...)
        elseif up == "RXX"
            # _RXX(ϕ) = exp(-i ϕ X⊗X) = exp(-i (2ϕ/2) X⊗X)
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), 2 * args[1]; kwargs...)
        elseif up == "RYY"
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), 2 * args[1]; kwargs...)
        elseif up == "RZZ"
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Z, t, :Z), args[1]; kwargs...)
        elseif up == "RXY"
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), args[1]; kwargs...)
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), args[1]; kwargs...)
        elseif up == "GIVENS"
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :X), args[1]; kwargs...)
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :Y), -args[1]; kwargs...)
        elseif up == "FSIM"
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), args[1]; kwargs...)
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), args[1]; kwargs...)
            return _camps_apply_cp!(st, q, t, args[2]; kwargs...)
        elseif up == "SWAPA"
            θ = π * args[1] / 2
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), θ; kwargs...)
            _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), θ; kwargs...)
            return _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Z, t, :Z), θ; kwargs...)
        end
    end

    # General 4×4 unitary fallback via exact KAK decomposition.
    A1, A2, θx, θy, θz, B1, B2 = _kak_decompose_2q(mat)
    _camps_apply_1q_unitary!(st, B1, q; kwargs...)
    _camps_apply_1q_unitary!(st, B2, t; kwargs...)
    abs(θx) > 1e-12 && _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :X, t, :X), θx; kwargs...)
    abs(θy) > 1e-12 && _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Y, t, :Y), θy; kwargs...)
    abs(θz) > 1e-12 && _camps_apply_pauli_string_rot!(st, _pauli_two(n, q, :Z, t, :Z), θz; kwargs...)
    _camps_apply_1q_unitary!(st, A1, q; kwargs...)
    _camps_apply_1q_unitary!(st, A2, t; kwargs...)
    return st
end

# Apply controlled 1-qubit unitary C-U(c, q) using fast paths + Barenco ZYZ decomposition.
function _camps_apply_controlled_1q!(st::CAMPSState, mat::AbstractMatrix, c::Int, q::Int; kwargs...)
    M = Matrix{ComplexF64}(mat)
    isunitary(M) || throw(ArgumentError("controlled single-qubit gate matrix must be unitary"))

    if isapprox(M, gate.I; atol=1e-12, rtol=1e-12)
        return st
    elseif isapprox(M, gate.X; atol=1e-12, rtol=1e-12)
        return _tableau_left_gate!(st.tableau, "CX", c, q)
    elseif isapprox(M, gate.Y; atol=1e-12, rtol=1e-12)
        return _tableau_left_gate!(st.tableau, "CY", c, q)
    elseif isapprox(M, gate.Z; atol=1e-12, rtol=1e-12)
        return _tableau_left_gate!(st.tableau, "CZ", c, q)
    elseif abs(M[1, 2]) < 1e-12 && abs(M[2, 1]) < 1e-12
        # Diagonal M = diag(e^{iφ0}, e^{iφ1}) => controlled phase P(φ0) on c + CP(φ1 - φ0) on (c, q)
        φ0, φ1 = angle(M[1, 1]), angle(M[2, 2])
        abs(φ0) > 1e-12 && _camps_apply_axis_rot!(st, :Z, c, φ0; kwargs...)
        abs(φ1 - φ0) > 1e-12 && _camps_apply_cp!(st, c, q, φ1 - φ0; kwargs...)
        return st
    end

    # Barenco Lemma 4.3: U = e^{iα} RZ(β) RY(γ) RZ(δ)
    α, β, γ, δ = zyz_decomposition(M)
    abs(α) > 1e-12 && _camps_apply_axis_rot!(st, :Z, c, α; kwargs...)
    abs(δ - β) > 1e-12 && _camps_apply_axis_rot!(st, :Z, q, (δ - β) / 2; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", c, q)
    abs(δ + β) > 1e-12 && _camps_apply_axis_rot!(st, :Z, q, -(δ + β) / 2; kwargs...)
    abs(γ) > 1e-12 && _camps_apply_axis_rot!(st, :Y, q, -γ / 2; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", c, q)
    abs(γ) > 1e-12 && _camps_apply_axis_rot!(st, :Y, q, γ / 2; kwargs...)
    abs(β) > 1e-12 && _camps_apply_axis_rot!(st, :Z, q, β; kwargs...)
    return st
end

function _camps_apply_ccz!(st::CAMPSState, a::Int, b::Int, c::Int; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", b, c)
    _camps_apply_axis_rot!(st, :Z, c, -π / 4; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", a, c)
    _camps_apply_axis_rot!(st, :Z, c, π / 4; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", b, c)
    _camps_apply_axis_rot!(st, :Z, c, -π / 4; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", a, c)
    _camps_apply_axis_rot!(st, :Z, b, π / 4; kwargs...)
    _camps_apply_axis_rot!(st, :Z, c, π / 4; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", a, b)
    _camps_apply_axis_rot!(st, :Z, a, π / 4; kwargs...)
    _camps_apply_axis_rot!(st, :Z, b, -π / 4; kwargs...)
    _tableau_left_gate!(st.tableau, "CX", a, b)
    return st
end

function _camps_apply_controlled_2q!(st::CAMPSState, mat::AbstractMatrix, c::Int, q::Int, t::Int; kwargs...)
    if isapprox(mat, gate.CZ; atol=1e-12, rtol=1e-12)
        return _camps_apply_ccz!(st, c, q, t; kwargs...)
    elseif isapprox(mat, gate.CX; atol=1e-12, rtol=1e-12)
        _tableau_left_gate!(st.tableau, "H", t)
        _camps_apply_ccz!(st, c, q, t; kwargs...)
        _tableau_left_gate!(st.tableau, "H", t)
        return st
    elseif isapprox(mat, gate.CY; atol=1e-12, rtol=1e-12)
        _tableau_left_gate!(st.tableau, "SD", t)
        _tableau_left_gate!(st.tableau, "H", t)
        _camps_apply_ccz!(st, c, q, t; kwargs...)
        _tableau_left_gate!(st.tableau, "H", t)
        _tableau_left_gate!(st.tableau, "S", t)
        return st
    elseif isapprox(mat, gate.SWAP; atol=1e-12, rtol=1e-12)
        _tableau_left_gate!(st.tableau, "CX", t, q)
        _tableau_left_gate!(st.tableau, "H", t)
        _camps_apply_ccz!(st, c, q, t; kwargs...)
        _tableau_left_gate!(st.tableau, "H", t)
        _tableau_left_gate!(st.tableau, "CX", t, q)
        return st
    else
        throw(ArgumentError("CAMPS supports controlled CX, CY, CZ, and SWAP as 3-qubit gates; decompose other 3-qubit operations into 1- and 2-qubit gates"))
    end
end

# =============================================================================
# Pauli Measurements, Conditional Branches, Reset & Noise Channels
# =============================================================================

# Projectively measure an arbitrary Hermitian Pauli string P on st = C|φ⟩.
# Mutates `st` in place and returns bit ∈ {0, 1} for eigenvalue (+1, -1).
function _camps_born_measure_pauli!(st::CAMPSState, p::_PauliString;
                                    maxdim::Union{Nothing,Int}=nothing,
                                    cutoff::Real=1e-12,
                                    disentangle::Bool=true,
                                    nsweeps::Int=2)
    p_tilde = _tableau_pullback(st.tableau, p)
    pivot = disentangle ? _camps_reduce_pauli_via_free_zeros!(st.tableau, st.mps, p_tilde) : 0

    # Gottesman-Knill fast path: P̃ collapsed to ±X_pivot on a free |0⟩ site!
    # Outcome is uniformly random, and absorbing H (or H*Z) keeps |φ⟩ unchanged.
    if pivot != 0
        bit = rand(Bool) ? 1 : 0
        eig = (1 - 2 * bit) * Int(p_tilde.phase)
        if eig == -1
            _tableau_right_dag_gate!(st.tableau, "Z", pivot)
        end
        _tableau_right_dag_gate!(st.tableau, "H", pivot)
        return bit
    end

    pauli_ops = _pauli_to_ops(p_tilde)
    if isempty(pauli_ops)
        # P̃ = ±I: deterministic outcome, state unchanged.
        return p_tilde.phase == 1 ? 0 : 1
    elseif length(pauli_ops) == 1
        op1 = only(pauli_ops)
        sites = it.siteinds(st.mps)
        q = op1.qubit
        # Rotate site q to Z basis, measure locally, rotate back.
        rot = op1.name == "X" ? gate.H : op1.name == "Y" ? gate.HSP : gate.I
        if op1.name != "Z"
            st.mps = la.normalize(it.apply(it.op(rot, sites[q]), st.mps; cutoff=cutoff))
        end
        st.mps, local_bit = born_measure_Z(st.mps, q)
        if op1.name != "Z"
            st.mps = la.normalize(it.apply(it.op(rot', sites[q]), st.mps; cutoff=cutoff))
        end
        # Convert local eigenvalue (+1/-1) to physical eigenvalue via p_tilde.phase.
        phys_bit = p_tilde.phase == 1 ? local_bit : (1 - local_bit)
        # If site q is now in an eigenstate of ±X/±Y/±Z, absorb a 1q Clifford to reset it to |0⟩!
        if disentangle
            if op1.name == "Z"
                local_bit == 1 && _camps_absorb_1q_clifford!(st, "X", q)
            elseif op1.name == "X"
                local_bit == 1 && _camps_absorb_1q_clifford!(st, "Z", q)
                _camps_absorb_1q_clifford!(st, "H", q)
            elseif op1.name == "Y"
                _camps_absorb_1q_clifford!(st, "SD", q)
                local_bit == 1 && _camps_absorb_1q_clifford!(st, "Z", q)
                _camps_absorb_1q_clifford!(st, "H", q)
            end
        end
        return phys_bit
    end

    p_phi = _camps_apply_pauli_ops_mps(st.mps, pauli_ops)
    expval = Int(p_tilde.phase) * real(it.inner(st.mps, p_phi))
    p0 = clamp((1.0 + expval) / 2.0, 0.0, 1.0)
    bit = rand() < p0 ? 0 : 1
    sign = (1 - 2 * bit) * Int(p_tilde.phase)
    projected = +(st.mps, sign * p_phi; alg="directsum")

    if disentangle
        pre_kw = _camps_trunc_kwargs(maxdim, cutoff; scale=2)
        it.truncate!(projected; pre_kw...)
        st.mps = la.normalize(projected)
        _camps_disentangle_sweep!(st; maxdim=maxdim, cutoff=cutoff, nsweeps=nsweeps)
    else
        trunc_kw = _camps_trunc_kwargs(maxdim, cutoff)
        it.truncate!(projected; trunc_kw...)
        st.mps = la.normalize(projected)
    end
    return bit
end

# Apply elementary 1q Clifford D to |φ⟩ and C ← C D† so C|φ⟩ is unchanged.
function _camps_absorb_1q_clifford!(st::CAMPSState, name::String, q::Int)
    sites = it.siteinds(st.mps)
    st.mps = la.normalize(it.apply(it.op(gates(name), sites[q]), st.mps; cutoff=0.0))
    _tableau_right_dag_gate!(st.tableau, name, q)
    return st
end

function _camps_born_measure!(st::CAMPSState, name::AbstractString, qubit::Int; kwargs...)
    rname = _resolve_measurement_name(String(name))
    axis = (rname == "MX" || rname == "M(X)") ? :X :
           (rname == "MY" || rname == "M(Y)") ? :Y : :Z
    return _camps_born_measure_pauli!(st, _pauli_single(st.tableau.n, qubit, axis); kwargs...)
end

function _born_measure(state::CAMPSState, o::QuantumOps; kwargs...)
    st = copy(state)
    bit = _camps_born_measure!(st, o.name, o.qubit; kwargs...)
    return st, bit
end

"""
    born_measure_Z(state::CAMPSState, qubit::Int;
                   cutoff=1e-12, maxdim=nothing, disentangle=true, nsweeps=2)

Sample a Z measurement of one qubit in a `CAMPSState`. Return the normalized
`(collapsed_state, bit)`, with bit 0 for `|0⟩` and 1 for `|1⟩`.
"""
function born_measure_Z(state::CAMPSState, qubit::Int;
                        cutoff::Real=1e-12,
                        maxdim::Union{Nothing,Int}=nothing,
                        disentangle::Bool=true,
                        nsweeps::Int=2)
    1 <= qubit <= get_N(state) || throw(ArgumentError("invalid qubit index $qubit"))
    st = copy(state)
    bit = _camps_born_measure_pauli!(st, _pauli_single(st.tableau.n, qubit, :Z);
                                     cutoff=cutoff, maxdim=maxdim,
                                     disentangle=disentangle, nsweeps=nsweeps)
    return st, bit
end

function __ifOp_apply(state::CAMPSState, name::String, qubit::Int,
                      if0::Vector{Op}, if1::Vector{Op}, noise=false;
                      track_measurements::Bool=false, kwargs...)
    camps_kw = _extract_camps_kwargs(kwargs)
    st = copy(state)
    bit = _camps_born_measure!(st, name, qubit; camps_kw...)
    return _apply_conditional_branch(st, bit, if0, if1, noise; track_measurements=track_measurements, camps_kw...)
end

# Pauli-basis expansion of a k-qubit matrix M (k ∈ {1, 2, 3}) acting on `qubits`.
const _PAULI_MATS_1Q = (gate.I, gate.X, gate.Z, gate.Y)
const _PAULI_SYMS_1Q = (:I, :X, :Z, :Y)

function _matrix_pauli_expansion(M::AbstractMatrix, n::Int, qubits::Tuple{Vararg{Int}}; tol::Real=1e-13)
    k = length(qubits)
    dim = 1 << k
    size(M) == (dim, dim) || throw(DimensionMismatch("expected $(dim)×$(dim) matrix for $k qubits"))
    terms = Tuple{ComplexF64,_PauliString}[]
    for idx in Iterators.product(ntuple(_ -> 1:4, k)...)
        Pmat = _PAULI_MATS_1Q[idx[1]]
        for j in 2:k
            Pmat = kron(Pmat, _PAULI_MATS_1Q[idx[j]])
        end
        coeff = la.tr(M * Pmat) / dim
        if abs(coeff) > tol
            x = falses(n)
            z = falses(n)
            for j in 1:k
                sym = _PAULI_SYMS_1Q[idx[j]]
                q = qubits[j]
                (sym === :X || sym === :Y) && (x[q] = true)
                (sym === :Z || sym === :Y) && (z[q] = true)
            end
            push!(terms, (ComplexF64(coeff), _PauliString(Int8(1), x, z)))
        end
    end
    return terms
end

function _camps_apply_general_kraus_op!(st::CAMPSState, K::AbstractMatrix, qubits::Tuple{Vararg{Int}};
                                        maxdim::Union{Nothing,Int}=nothing,
                                        cutoff::Real=1e-12,
                                        disentangle::Bool=true,
                                        nsweeps::Int=2)
    # If K is a 2-qubit product K = A ⊗ B, apply A then B sequentially.
    if length(qubits) == 2
        C = reshape(permutedims(reshape(Matrix{ComplexF64}(K), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
        F = la.svd(C)
        if F.S[1] > 1e-14 && (length(F.S) == 1 || F.S[2] <= 1e-12 * F.S[1])
            A = sqrt(F.S[1]) * reshape(F.Vt[1, :], 2, 2)
            B = sqrt(F.S[1]) * reshape(F.U[:, 1], 2, 2)
            _camps_apply_general_kraus_op!(st, A, (qubits[1],); maxdim=maxdim, cutoff=cutoff, disentangle=disentangle, nsweeps=nsweeps)
            _camps_apply_general_kraus_op!(st, B, (qubits[2],); maxdim=maxdim, cutoff=cutoff, disentangle=disentangle, nsweeps=nsweeps)
            return st
        end
    end

    terms = _matrix_pauli_expansion(K, st.tableau.n, qubits)
    isempty(terms) && throw(ArgumentError("selected Kraus operator is identically zero"))

    # Sum c_P * (C† P C) |φ⟩ via MPS direct sum.
    acc = nothing
    for (coeff, p) in terms
        p_tilde = _tableau_pullback(st.tableau, p)
        weight = coeff * Int(p_tilde.phase)
        term_mps = weight * _camps_apply_pauli_ops_mps(st.mps, _pauli_to_ops(p_tilde))
        acc = acc === nothing ? term_mps : +(acc, term_mps; alg="directsum")
    end

    if disentangle
        pre_kw = _camps_trunc_kwargs(maxdim, cutoff; scale=2)
        it.truncate!(acc; pre_kw...)
        st.mps = la.normalize(acc)
        _camps_disentangle_sweep!(st; maxdim=maxdim, cutoff=cutoff, nsweeps=nsweeps)
    else
        trunc_kw = _camps_trunc_kwargs(maxdim, cutoff)
        it.truncate!(acc; trunc_kw...)
        st.mps = la.normalize(acc)
    end
    return st
end

function _is_reset_kraus(kraus_ops::Vector{<:AbstractMatrix})
    length(kraus_ops) == 2 || return false
    size(kraus_ops[1]) == (2, 2) && size(kraus_ops[2]) == (2, 2) || return false
    return isapprox(kraus_ops[1], ComplexF64[1 0; 0 0]; atol=1e-12, rtol=0) &&
           isapprox(kraus_ops[2], ComplexF64[0 1; 0 0]; atol=1e-12, rtol=0)
end

function _camps_apply_quantum_channel!(st::CAMPSState, kraus_ops::Vector{<:AbstractMatrix},
                                       qubits::Int...; kwargs...)
    # Fast path for single-qubit ideal reset |0⟩⟨0|, |0⟩⟨1|.
    if length(qubits) == 1 && _is_reset_kraus(kraus_ops)
        bit = _camps_born_measure_pauli!(st, _pauli_single(st.tableau.n, qubits[1], :Z); kwargs...)
        bit == 1 && _tableau_left_gate!(st.tableau, "X", qubits[1])
        return st
    end

    # Fast path for unitary-proportional Kraus channels (Pauli & coherent rotation noise).
    grams = [K' * K for K in kraus_ops]
    weights = [real(la.tr(G)) / size(G, 1) for G in grams]
    tolerances = [8eps(real(float(one(eltype(G))))) for G in grams]
    if all(isapprox(G, p * la.I; rtol=tol, atol=0) for (G, p, tol) in zip(grams, weights, tolerances))
        idx = _weighted_sample(weights)
        w = weights[idx]
        U = kraus_ops[idx] / sqrt(w)
        if length(qubits) == 1
            return _camps_apply_1q_unitary!(st, U, qubits[1]; kwargs...)
        elseif length(qubits) == 2
            return _camps_apply_2q_unitary!(st, U, qubits[1], qubits[2]; kwargs...)
        else
            return _camps_apply_general_kraus_op!(st, U, qubits; kwargs...)
        end
    end

    # General non-unitary Kraus channel: compute exact Born probabilities from K†K.
    probs = Vector{Float64}(undef, length(kraus_ops))
    for (m, G) in enumerate(grams)
        p_m = 0.0
        for (coeff, P) in _matrix_pauli_expansion(G, st.tableau.n, qubits)
            p_m += real(coeff) * _camps_expect_pauli(st, P)
        end
        probs[m] = max(0.0, p_m)
    end
    idx = _weighted_sample(probs)
    return _camps_apply_general_kraus_op!(st, kraus_ops[idx], qubits; kwargs...)
end

function apply_noise(state::CAMPSState, op::QuantumOps, noise::NoiseModel; kwargs...)
    (hasproperty(op, :noisy) && op.noisy) || return state
    st = copy(state)
    _camps_apply_noise!(st, op, noise; kwargs...)
    return st
end

function _camps_apply_noise!(st::CAMPSState, op::QuantumOps, noise::NoiseModel; kwargs...)
    (hasproperty(op, :noisy) && op.noisy) || return st
    if op.q == 1
        if op.control == -2
            if noise.q1 isa QuantumChannel
                _camps_apply_quantum_channel!(st, noise.q1.kraus, op.qubit; kwargs...)
            end
        elseif noise.q2 isa QuantumChannel
            _camps_apply_quantum_channel!(st, noise.q2.kraus, op.control, op.qubit; kwargs...)
        end
    elseif op.q == 2 && noise.q2 isa QuantumChannel
        _camps_apply_quantum_channel!(st, noise.q2.kraus, op.qubit, op.target_qubit; kwargs...)
    end
    return st
end

# =============================================================================
# Core apply & camps Entry Points
# =============================================================================

function _extract_camps_kwargs(kwargs)
    cutoff = get(kwargs, :cutoff, 1e-12)
    maxdim = get(kwargs, :maxdim, nothing)
    disentangle = get(kwargs, :disentangle, true)
    nsweeps = get(kwargs, :nsweeps, 2)
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    (maxdim === nothing || (maxdim isa Integer && maxdim > 0)) ||
        throw(ArgumentError("maxdim must be nothing or a positive integer"))
    nsweeps >= 0 || throw(ArgumentError("nsweeps must be nonnegative"))
    return (cutoff=cutoff, maxdim=maxdim === nothing ? nothing : Int(maxdim),
            disentangle=Bool(disentangle), nsweeps=Int(nsweeps))
end

function _camps_apply_op!(st::CAMPSState, op::QuantumOps;
                          noise::Union{NoiseModel,Bool}=false,
                          track_measurements::Bool=false,
                          pars=nothing,
                          cutoff::Real=1e-12,
                          maxdim::Union{Nothing,Int}=nothing,
                          disentangle::Bool=true,
                          nsweeps::Int=2)
    mid_measurements = Int[]
    ckw = (cutoff=cutoff, maxdim=maxdim, disentangle=disentangle, nsweeps=nsweeps)

    if op isa OpF
        if op.data isa Vector{<:QuantumOps}
            for nested in op.data
                st, nested_bits = _camps_apply_op!(st, nested; noise=noise,
                                                   track_measurements=track_measurements,
                                                   pars=pars, ckw...)
                track_measurements && append!(mid_measurements, nested_bits)
            end
        elseif op.data isa Function
            st = op.data(st)
            st isa CAMPSState || throw(ArgumentError("OpF callback must return a CAMPSState"))
        elseif op.data isa AbstractMatrixS
            dim = size(op.data, 1)
            if dim == 2 && get_N(st) == 1
                _camps_apply_1q_unitary!(st, op.data, 1; ckw...)
            elseif dim == 4 && get_N(st) == 2
                _camps_apply_2q_unitary!(st, op.data, 1, 2; ckw...)
            else
                throw(ArgumentError("dense full-register OpF matrices above 2 qubits are not supported on CAMPSState"))
            end
        end
    elseif op isa OpQC
        if op.q == 1
            _camps_apply_quantum_channel!(st, op.kraus, op.qubit; ckw...)
        elseif op.q == 2
            _camps_apply_quantum_channel!(st, op.kraus, op.qubit, op.target_qubit; ckw...)
        elseif op.q == 3
            _camps_apply_quantum_channel!(st, op.kraus, op.qubit, op.qubit + 1, op.qubit + 2; ckw...)
        else
            throw(ArgumentError("unsupported OpQC arity $(op.q) on CAMPSState"))
        end
    elseif op.type == "🔬"
        if op isa ifOp
            st, record = op.born_apply(st, noise; track_measurements=track_measurements, ckw...)
            track_measurements && append!(mid_measurements, record)
        else
            bit = _camps_born_measure!(st, op.name, op.qubit; ckw...)
            track_measurements && push!(mid_measurements, bit)
        end
    else
        mat = op.type == "f" ? op.mat(pars...) : op.mat
        if op.control == -2
            if op.q == 1
                up = uppercase(_clean_name(op.name))
                if up in ("I", "X", "Y", "Z", "H", "S", "SD", "SDG", "SX", "XSQRT",
                          "SXDG", "SXD", "XSQRTD", "XSQRTDG", "HSP", "HY") &&
                   isapprox(mat, gates(up); atol=1e-12, rtol=1e-12)
                    _tableau_left_gate!(st.tableau, up, op.qubit)
                elseif up == "T" && isapprox(mat, gate.T; atol=1e-12, rtol=1e-12)
                    _camps_apply_axis_rot!(st, :Z, op.qubit, π / 4; ckw...)
                elseif (up == "TD" || up == "TDG") && isapprox(mat, gate.TD; atol=1e-12, rtol=1e-12)
                    _camps_apply_axis_rot!(st, :Z, op.qubit, -π / 4; ckw...)
                else
                    _camps_apply_1q_unitary!(st, mat, op.qubit; ckw...)
                end
            elseif op.q == 2
                _camps_apply_2q_unitary!(st, mat, op.qubit, op.target_qubit; name=op.name, ckw...)
            else
                throw(ArgumentError("unsupported $(op.q)-qubit gate $(op.name) on CAMPSState"))
            end
        else
            if op.q == 1
                _camps_apply_controlled_1q!(st, mat, op.control, op.qubit; ckw...)
            elseif op.q == 2
                _camps_apply_controlled_2q!(st, mat, op.control, op.qubit, op.target_qubit; ckw...)
            else
                throw(ArgumentError("unsupported controlled $(op.q)-qubit gate $(op.name) on CAMPSState"))
            end
        end
    end

    if noise isa NoiseModel
        _camps_apply_noise!(st, op, noise; ckw...)
    end

    return st, mid_measurements
end

function apply(
    state::CAMPSState,
    op::QuantumOps;
    noise::Union{NoiseModel,Bool}=false,
    track_measurements::Bool=false,
    kwargs...
)
    _unsupported_keyword_check(kwargs, [:cutoff, :maxdim, :disentangle, :nsweeps, :pars])
    ckw = _extract_camps_kwargs(kwargs)
    pars = get(kwargs, :pars, nothing)
    st = copy(state)
    st, mid_measurements = _camps_apply_op!(st, op; noise=noise,
                                            track_measurements=track_measurements,
                                            pars=pars, ckw...)
    return track_measurements ? (st, mid_measurements) : st
end

function apply(
    state::CAMPSState,
    vector_of_ops::Union{Vector{<:QuantumOps},AbstractVector{<:Tuple}};
    noise::Union{NoiseModel,Bool}=false,
    track_measurements::Bool=false,
    kwargs...
)
    _unsupported_keyword_check(kwargs, [:cutoff, :maxdim, :disentangle, :nsweeps, :pars])
    ckw = _extract_camps_kwargs(kwargs)
    pars = get(kwargs, :pars, nothing)
    st = copy(state)
    mid_measurements = Int[]

    if vector_of_ops isa Vector{<:QuantumOps}
        for op in vector_of_ops
            st, new_vals = _camps_apply_op!(st, op; noise=noise,
                                            track_measurements=track_measurements,
                                            pars=pars, ckw...)
            track_measurements && append!(mid_measurements, new_vals)
        end
    else
        for tup in vector_of_ops
            st, new_vals = _camps_apply_op!(st, Op(tup...); noise=noise,
                                            track_measurements=track_measurements,
                                            pars=pars, ckw...)
            track_measurements && append!(mid_measurements, new_vals)
        end
    end

    return track_measurements ? (st, mid_measurements) : st
end

"""
    camps(ops; N=get_N(ops), noise=false, cutoff=1e-12, maxdim=nothing,
          disentangle=true, nsweeps=2) -> CAMPSState

Evolve the `N`-qubit zero state `|0…0⟩` through `ops` using the
Clifford-Augmented Matrix Product State (`CAMPSState`) backend.

Clifford gates are absorbed into the `CliffordTableau` in `O(N)` time with zero
MPS bond-dimension growth. Non-Clifford gates pull their Pauli generators
through the Clifford frame, apply free-zero reduction and two-qubit Clifford
disentangling sweeps (`disentangle=true`, `nsweeps=2`), and truncate the MPS
core using `cutoff` and optional `maxdim`.
"""
function camps(ops::Union{Vector{<:QuantumOps},AbstractVector{<:Tuple}};
               N::Int=(ops isa Vector{<:QuantumOps} ? get_N(ops) : get_N(Op[Op(t...) for t in ops])),
               noise::Union{NoiseModel,Bool}=false,
               kwargs...)
    N > 0 || throw(ArgumentError("N must be positive; supply N explicitly for an empty operation list"))
    ops_n = ops isa Vector{<:QuantumOps} ? get_N(ops) : get_N(Op[Op(t...) for t in ops])
    N >= ops_n || throw(ArgumentError("N must include every operated qubit"))
    return apply(CAMPSState(N), ops; noise=noise, kwargs...)
end

# =============================================================================
# Diagnostics, Expectations, Sampling, Conversion & Fidelity
# =============================================================================

la.norm(st::CAMPSState) = la.norm(st.mps)
la.normalize(st::CAMPSState) = CAMPSState(copy(st.tableau), la.normalize(st.mps))

maxlinkdim(st::CAMPSState) = length(st.mps) <= 1 ? 1 : it.maxlinkdim(st.mps)

entanglement_entropy(st::CAMPSState, b::Int; spectrum_bool::Bool=false) =
    entanglement_entropy(st.mps, b; spectrum_bool=spectrum_bool)
entanglement_entropy(st::CAMPSState; spectrum_bool::Bool=false) =
    entanglement_entropy(st.mps; spectrum_bool=spectrum_bool)

function _camps_expect_pauli(st::CAMPSState, p::_PauliString)
    p_tilde = _tableau_pullback(st.tableau, p)
    pauli_ops = _pauli_to_ops(p_tilde)
    isempty(pauli_ops) && return Float64(p_tilde.phase)
    val = real(it.inner(st.mps, _camps_apply_pauli_ops_mps(st.mps, pauli_ops)))
    return Float64(p_tilde.phase) * val
end

function expect(st::CAMPSState, op_str::String, qubit::Int)
    n = get_N(st)
    1 <= qubit <= n || throw(ArgumentError("qubit index $qubit out of range 1:$n"))
    up = uppercase(strip(op_str))
    if up == "I"
        return 1.0
    elseif up == "X"
        return _camps_expect_pauli(st, _pauli_single(n, qubit, :X))
    elseif up == "Y"
        return _camps_expect_pauli(st, _pauli_single(n, qubit, :Y))
    elseif up == "Z"
        return _camps_expect_pauli(st, _pauli_single(n, qubit, :Z))
    end
    mat = gates(op_str)
    size(mat) == (2, 2) || throw(ArgumentError("expected a single-qubit operator, got `$op_str`"))
    total = 0.0
    for (coeff, P) in _matrix_pauli_expansion(mat, n, (qubit,))
        total += real(coeff) * _camps_expect_pauli(st, P)
    end
    return total
end

function expect(st::CAMPSState, op_str::String)
    n = get_N(st)
    if occursin(',', op_str)
        parts = strip.(split(op_str, ','))
        length(parts) == n || throw(DimensionMismatch("expected $n comma-separated operators, got $(length(parts))"))
        return correlation(st, op_str, collect(1:n))
    end
    return [expect(st, op_str, q) for q in 1:n]
end

function expect(st::CAMPSState, op::Op)
    n = get_N(st)
    op.mat isa AbstractMatrix || throw(ArgumentError("expect requires a bound matrix Op"))
    if op.control == -2 && op.q == 1
        total = 0.0
        for (coeff, P) in _matrix_pauli_expansion(op.mat, n, (op.qubit,))
            total += real(coeff) * _camps_expect_pauli(st, P)
        end
        return total
    elseif op.control == -2 && op.q == 2
        total = 0.0
        for (coeff, P) in _matrix_pauli_expansion(op.mat, n, (op.qubit, op.target_qubit))
            total += real(coeff) * _camps_expect_pauli(st, P)
        end
        return total
    elseif op.control > 0 && op.q == 1
        full = kron(gate.P0, Matrix{ComplexF64}(la.I, 2, 2)) + kron(gate.P1, Matrix{ComplexF64}(op.mat))
        total = 0.0
        for (coeff, P) in _matrix_pauli_expansion(full, n, (op.control, op.qubit))
            total += real(coeff) * _camps_expect_pauli(st, P)
        end
        return total
    elseif op.control > 0 && op.q == 2
        full = kron(gate.P0, Matrix{ComplexF64}(la.I, 4, 4)) + kron(gate.P1, Matrix{ComplexF64}(op.mat))
        total = 0.0
        for (coeff, P) in _matrix_pauli_expansion(full, n, (op.control, op.qubit, op.target_qubit))
            total += real(coeff) * _camps_expect_pauli(st, P)
        end
        return total
    end
    throw(ArgumentError("unsupported Op in expect(::CAMPSState, ::Op)"))
end

function expect(st::CAMPSState, matrix::AbstractMatrixS)
    vec = to_state(st)
    return real(vec' * matrix * vec)
end

function correlation(st::CAMPSState, names::String, qubits::Vector)
    n = get_N(st)
    operators = strip.(split(names, ','))
    length(operators) == length(qubits) || throw(DimensionMismatch("one operator is required per qubit"))
    all(q -> q isa Integer && 1 <= q <= n, qubits) || throw(ArgumentError("qubit indices must lie in 1:$n"))

    # Fast path when all factors are Paulis on distinct qubits.
    ups = uppercase.(String.(operators))
    if allunique(qubits) && all(u -> u in ("I", "X", "Y", "Z"), ups)
        x = falses(n)
        z = falses(n)
        for (u, q) in zip(ups, qubits)
            (u == "X" || u == "Y") && (x[q] = true)
            (u == "Z" || u == "Y") && (z[q] = true)
        end
        return _camps_expect_pauli(st, _PauliString(Int8(1), x, z))
    end

    # General product of 1-qubit operators (including projectors or repeated qubits).
    terms = Tuple{ComplexF64,_PauliString}[(1.0 + 0.0im, _PauliString(Int8(1), falses(n), falses(n)))]
    for (op_name, q) in zip(operators, qubits)
        mat = gates(String(op_name))
        size(mat) == (2, 2) || throw(ArgumentError("correlation expects single-qubit operators, got `$op_name`"))
        local_terms = _matrix_pauli_expansion(mat, n, (Int(q),))
        next_terms = Tuple{ComplexF64,_PauliString}[]
        for (c1, p1) in terms, (c2, p2) in local_terms
            ipow, p_prod = _pauli_mul_raw(p1, p2)
            c_prod = c1 * c2 * (ipow == 1 ? im : 1.0)
            abs(c_prod) > 1e-14 && push!(next_terms, (c_prod, p_prod))
        end
        terms = next_terms
    end

    total = 0.0 + 0.0im
    for (c, p) in terms
        total += c * _camps_expect_pauli(st, p)
    end
    return real(total)
end

correlation(st::CAMPSState, qubits::Vector) = correlation(st, join(fill("Z", length(qubits)), ','), qubits)

# Fast stabilizer Z-basis shot when core MPS is |0…0⟩ (zero ITensor calls per shot).
function _camps_sample_stabilizer_shot(base_pz::Vector{_PauliString})
    n = length(base_pz)
    pz = [copy(p) for p in base_pz]
    bits = Vector{Int}(undef, n)

    for j in 1:n
        p = pz[j]
        pivot = findfirst(p.x)
        if pivot === nothing
            bits[j] = p.phase == 1 ? 0 : 1
        else
            # Reduce p to ±X_pivot via right-conjugations on remaining pz[j+1:n].
            if p.z[pivot]
                _pauli_conj_gate!(p, "SD", pivot)
                for k in (j + 1):n
                    _pauli_conj_gate!(pz[k], "SD", pivot)
                end
            end
            for a in 1:n
                a == pivot && continue
                xa, za = p.x[a], p.z[a]
                if xa || za
                    gname = (xa && za) ? "CY" : (xa ? "CX" : "CZ")
                    _pauli_conj_gate!(p, gname, pivot, a)
                    for k in (j + 1):n
                        _pauli_conj_gate!(pz[k], gname, pivot, a)
                    end
                end
            end
            b = rand(Bool) ? 1 : 0
            bits[j] = b
            eig = (1 - 2 * b) * Int(p.phase)
            if eig == -1
                for k in (j + 1):n
                    _pauli_conj_gate!(pz[k], "Z", pivot)
                end
            end
            for k in (j + 1):n
                _pauli_conj_gate!(pz[k], "H", pivot)
            end
        end
    end
    return bits
end

function sample_bit(state::CAMPSState, shots::Int=1; basis::AbstractString="Z", cutoff::Real=1e-12)
    shots >= 0 || throw(ArgumentError("shots must be nonnegative"))
    shots == 0 && return Vector{Int}[]
    n = get_N(state)
    up_basis = uppercase(String(basis))
    st = copy(state)
    if up_basis == "X"
        for q in 1:n
            _tableau_left_gate!(st.tableau, "H", q)
        end
    elseif up_basis != "Z"
        throw(ArgumentError("basis must be Z or X"))
    end

    # Fast stabilizer path when the MPS core is the product state |0…0⟩.
    z_exps = real.(it.expect(st.mps, "Z"))
    if all(z -> z > 1.0 - 1e-11, z_exps)
        return [_camps_sample_stabilizer_shot(st.tableau.pz) for _ in 1:shots]
    end

    results = Vector{Vector{Int}}(undef, shots)
    for s in 1:shots
        shot = copy(st)
        bits = Vector{Int}(undef, n)
        for j in 1:n
            bits[j] = _camps_born_measure_pauli!(shot, _pauli_single(n, j, :Z);
                                                 cutoff=cutoff, disentangle=true, nsweeps=1)
        end
        results[s] = bits
    end
    return results
end

sample(state::CAMPSState, shots::Int=1) = bin2int.(sample_bit(state, shots))

function sample_exact(state::CAMPSState; max_qubits::Int=16)
    n = get_N(state)
    n <= max_qubits ||
        throw(ArgumentError("exact CAMPS probabilities require 2^N outcomes; use measure(state, shots), or explicitly raise max_qubits"))
    return collect(0:(2^n - 1)), Vector(abs2.(to_state(state)))
end

function prob(state::CAMPSState, shots::Int, select_qubit::Int)
    shots > 0 || throw(ArgumentError("shots must be positive"))
    1 <= select_qubit <= get_N(state) || throw(ArgumentError("invalid selected qubit"))
    all_sample = hcat(sample_bit(state, shots)...)[select_qubit, :]
    P1 = sum(all_sample) / shots
    return 1 - P1, P1
end

"""
    to_state(state::CAMPSState)

Materialize `|Ψ⟩ = C |φ⟩` as a sparse `2^N` statevector (with qubit 1 the most
significant bit), up to a single global phase. Contracts the core MPS first and
then applies the synthesized Clifford circuit `clifford_ops(state)` directly to
the vector.
"""
function to_state(state::CAMPSState)
    vec = to_state(state.mps)
    _tableau_is_identity(state.tableau) && return vec
    return apply(vec, clifford_ops(state.tableau))
end

"""
    to_rho(state::CAMPSState)

Materialize the exact `2^N × 2^N` pure-state density matrix `|Ψ⟩⟨Ψ|`. Because
global phases cancel in `|Ψ⟩⟨Ψ|`, `to_rho(state)` is independent of the global
phase of `clifford_ops(state)`.
"""
function to_rho(state::CAMPSState)
    vec = to_state(state)
    return vec * vec'
end

"""
    to_MPS(state::CAMPSState; kwargs...)

Convert a `CAMPSState` into a standard `ITensorMPS.MPS` representing `C |φ⟩`
(up to a single global phase) by applying `clifford_ops(state)` to `state.mps`.
Pass `cutoff` or `maxdim` via `kwargs` to control MPS truncation during
Clifford re-absorption; note that re-absorbing a deeply entangling Clifford
frame can substantially increase the MPS bond dimension.
"""
function to_MPS(state::CAMPSState; kwargs...)
    psi = copy(state.mps)
    _tableau_is_identity(state.tableau) && return psi
    return apply(psi, clifford_ops(state.tableau); kwargs...)
end

function fidelity(a::CAMPSState, b::CAMPSState)
    get_N(a) == get_N(b) || throw(DimensionMismatch("states must have the same number of qubits"))
    if a.tableau == b.tableau
        return fidelity(a.mps, b.mps)
    elseif get_N(a) <= 16
        return fidelity(to_state(a), to_state(b))
    else
        rel = _tableau_relative(a.tableau, b.tableau)
        moved = apply(copy(b.mps), clifford_ops(rel); cutoff=1e-12)
        return fidelity(a.mps, moved)
    end
end

fidelity(a::CAMPSState, b::AbstractVectorS) = fidelity(to_state(a), b)
fidelity(a::AbstractVectorS, b::CAMPSState) = fidelity(a, to_state(b))

function fidelity(a::CAMPSState, b::it.MPS)
    get_N(a) == get_N(b) || throw(DimensionMismatch("states must have the same number of qubits"))
    return get_N(a) <= 16 ? fidelity(to_state(a), to_state(b)) : fidelity(to_MPS(a; cutoff=1e-12), b)
end

fidelity(a::it.MPS, b::CAMPSState) = fidelity(b, a)

# =============================================================================
# QEC Integration (StabilizerCode & qec_round)
# =============================================================================

function _qec_append_ancillas(st::CAMPSState, m::Int)
    m == 0 && return copy(st)
    m > 0 || throw(ArgumentError("number of ancillas must be nonnegative"))
    n = get_N(st)
    N = n + m
    pad = falses(m)
    px = Vector{_PauliString}(undef, N)
    pz = Vector{_PauliString}(undef, N)
    for j in 1:n
        px[j] = _PauliString(st.tableau.px[j].phase, vcat(st.tableau.px[j].x, pad), vcat(st.tableau.px[j].z, pad))
        pz[j] = _PauliString(st.tableau.pz[j].phase, vcat(st.tableau.pz[j].x, pad), vcat(st.tableau.pz[j].z, pad))
    end
    for i in 1:m
        px[n + i] = _pauli_single(N, n + i, :X)
        pz[n + i] = _pauli_single(N, n + i, :Z)
    end
    mps_ext = _qec_append_ancillas(st.mps, m)
    return CAMPSState(CliffordTableau(N, px, pz), mps_ext)
end

function _permute_pauli(p::_PauliString, perm::Vector{Int})
    n = length(p.x)
    x = falses(n)
    z = falses(n)
    @inbounds for i in 1:n
        x[perm[i]] = p.x[i]
        z[perm[i]] = p.z[i]
    end
    return _PauliString(p.phase, x, z)
end

function _qec_embed_logical_camps(logical_state::CAMPSState, n::Int, logical_indices::AbstractVector)
    k = get_N(logical_state)
    ext = _qec_append_ancillas(logical_state, n - k)
    ancillas = setdiff(1:n, logical_indices)
    # perm[old_site] = new_site: old 1:k -> logical_indices, old (k+1):n -> ancillas.
    perm = vcat(Int.(logical_indices), ancillas)
    px_new = Vector{_PauliString}(undef, n)
    pz_new = Vector{_PauliString}(undef, n)
    for old_i in 1:n
        px_new[perm[old_i]] = _permute_pauli(ext.tableau.px[old_i], perm)
        pz_new[perm[old_i]] = _permute_pauli(ext.tableau.pz[old_i], perm)
    end
    mps_new = it.movesites(ext.mps, collect(1:k) .=> Int.(logical_indices); cutoff=1e-14)
    return CAMPSState(CliffordTableau(n, px_new, pz_new), mps_new)
end

function reduce_and_relabel_qubits(state::CAMPSState, logical_indices::Vector, qubit_mapping::Vector=Int[])
    n = get_N(state)
    k = length(logical_indices)
    if !isempty(qubit_mapping)
        mapping = [qubit_mapping; (length(qubit_mapping) + 1):n]
        sort(mapping) == collect(1:n) || throw(ArgumentError("qubit_mapping must permute its qubits without repetition"))
        logical_indices = [findfirst(==(q), mapping) for q in logical_indices]
    end
    all(q -> q isa Integer && 1 <= q <= n, logical_indices) && allunique(logical_indices) ||
        throw(ArgumentError("retained qubits must be distinct valid indices"))

    discarded = setdiff(1:n, logical_indices)
    for q in discarded
        isapprox(expect(state, "Z", q), 1.0; atol=1e-10) ||
            throw(ArgumentError("discarded qubits are not in |0⟩; correct/reset them first"))
    end

    st = copy(state)
    sites = it.siteinds(st.mps)

    # Decouple every discarded qubit in both tableau and MPS so px[q]=+X_q, pz[q]=+Z_q,
    # and site q of st.mps is the product state |0⟩_q.
    locked = Set{Int}()
    for q in discarded
        function apply_both!(name::String, a::Int, b::Int=-1)
            op_t = b < 0 ? it.op(gates(name), sites[a]) : it.op(gates(name), sites[b], sites[a])
            st.mps = la.normalize(it.apply(op_t, st.mps; cutoff=1e-12))
            _tableau_right_dag_gate!(st.tableau, name, a, b)
        end

        pz = st.tableau.pz[q]
        support = [s for s in 1:n if !(s in locked) && (pz.x[s] || pz.z[s])]
        for s in support
            xs, zs = pz.x[s], pz.z[s]
            if xs && zs
                apply_both!("S", s)
                apply_both!("H", s)
            elseif xs
                apply_both!("H", s)
            end
        end
        pivot = first(support)
        for idx in 2:length(support)
            apply_both!("CX", support[idx], pivot)
        end
        if pivot != q
            apply_both!("SWAP", q, pivot)
        end
        if st.tableau.pz[q].phase == -1
            apply_both!("X", q)
        end

        # Now pz[q] == +Z_q => site q of st.mps is |0⟩_q!
        # Controlled Cliffords on q and SD(q)/Z(q) leave |0⟩_q unchanged, so update tableau only.
        px = st.tableau.px[q]
        if px.x[q] && px.z[q]
            _tableau_right_dag_gate!(st.tableau, "SD", q)
        end
        for a in 1:n
            a == q && continue
            xa, za = px.x[a], px.z[a]
            if xa && za
                _tableau_right_dag_gate!(st.tableau, "CY", q, a)
            elseif xa
                _tableau_right_dag_gate!(st.tableau, "CX", q, a)
            elseif za
                _tableau_right_dag_gate!(st.tableau, "CZ", q, a)
            end
        end
        if st.tableau.px[q].phase == -1
            _tableau_right_dag_gate!(st.tableau, "Z", q)
        end
        push!(locked, q)
    end

    # Extract retained qubits in the requested order.
    mps_reduced = reduce_and_relabel_qubits(st.mps, Int.(logical_indices))
    k == 0 && return mps_reduced

    px_red = Vector{_PauliString}(undef, k)
    pz_red = Vector{_PauliString}(undef, k)
    for (new_j, old_j) in enumerate(logical_indices)
        px_old = st.tableau.px[old_j]
        pz_old = st.tableau.pz[old_j]
        x_px = BitVector([px_old.x[q] for q in logical_indices])
        z_px = BitVector([px_old.z[q] for q in logical_indices])
        x_pz = BitVector([pz_old.x[q] for q in logical_indices])
        z_pz = BitVector([pz_old.z[q] for q in logical_indices])
        px_red[new_j] = _PauliString(px_old.phase, x_px, z_px)
        pz_red[new_j] = _PauliString(pz_old.phase, x_pz, z_pz)
    end
    return CAMPSState(CliffordTableau(k, px_red, pz_red), mps_reduced)
end

"""
    qec_round(code::StabilizerCode, state::CAMPSState;
              decoder=nothing, cutoff=1e-12, maxdim=nothing,
              disentangle=true, nsweeps=2)

Measure the stabilizers of `code` ideally on a `CAMPSState`, apply Pauli
recovery, and return `(state=corrected_state, syndrome=bits)`.

Stabilizer projectors pull directly through the Clifford tableau; on stabilizer
codewords (with arbitrary Pauli errors), the syndrome round completes in
`O(m N²)` tableau steps while keeping `maxlinkdim(state) == 1`. On non-Clifford
encoded states, `cutoff`, `maxdim`, `disentangle`, and `nsweeps` control MPS
compression and Clifford disentangling after each stabilizer projection.
"""
function qec_round(code::StabilizerCode, state::CAMPSState;
                   decoder=nothing,
                   cutoff::Real=1e-12,
                   maxdim::Union{Nothing,Int}=nothing,
                   disentangle::Bool=true,
                   nsweeps::Int=2)
    get_N(state) == code.n || throw(DimensionMismatch("an ideal QEC round requires $(code.n) data qubits"))
    state_norm = la.norm(state)
    isfinite(state_norm) && state_norm > 0 || throw(ArgumentError("the input state must have a finite, positive norm"))
    st = la.normalize(state)
    bits = Int[]
    n = code.n
    for row in eachrow(code.generator)
        x = BitVector(row[1:n] .== 1)
        z = BitVector(row[(n + 1):(2n)] .== 1)
        stab = _PauliString(Int8(1), x, z)
        bit = _camps_born_measure_pauli!(st, stab; cutoff=cutoff, maxdim=maxdim,
                                         disentangle=disentangle, nsweeps=nsweeps)
        push!(bits, bit)
    end
    return (state=code.correct(st, bits; decoder=decoder), syndrome=bits)
end
