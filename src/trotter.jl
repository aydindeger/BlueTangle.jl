"""
    evolve(state, H, time)

Exact evolution by `exp(-im*time*H)` for a statevector or density matrix.
This materializes the matrix exponential; use `hamiltonian_exp` for a circuit
approximation when the full Hilbert-space matrix is too large.
"""
function evolve(state::Union{AbstractVectorS,AbstractMatrixS}, H::AbstractMatrix, time::Real)
    size(H,1) == size(H,2) == size(state,1) || throw(DimensionMismatch("Hamiltonian and state dimensions must agree"))
    U = exp(-im * time * Matrix(H))
    result = state isa AbstractVector ? U * state : U * state * U'
    return sa.issparse(state) ? sa.sparse(result) : result
end

"""
    trotter_ising(N, total_time, J, h; dt=0.1, order=1)

Product-formula evolution for `H = -J*sum(Z_i*Z_{i+1}) - h*sum(X_i)`.
"""
trotter_ising(N::Int, total_time::Real, J::Real, h::Real; dt=0.1, order=1) =
    hamiltonian_exp(N, total_time, [-h, "X", -J, "Z,Z"]; dt=dt, order=order)

# exp(-i*angle/2*P): rotate each nonidentity Pauli to Z, collect its parity,
# rotate the parity qubit, and undo the parity circuit and basis rotations.
function _apply_term(all_ops::Vector{QuantumOps}, terms::AbstractVector, angle::Real, qubits::Vector{Int})
    iszero(angle) && return
    active = findall(!=("I"), terms)
    if isempty(active)
        push!(all_ops, Op("global phase", exp(-im*angle/2) * gate.I, first(qubits)))
        return
    end
    if length(active) == 1
        index = only(active)
        pauli = terms[index]
        matrix = pauli == "X" ? _RX(angle) : pauli == "Y" ? _RY(angle) : _RZ(angle)
        push!(all_ops, Op("R$pauli($angle)", matrix, qubits[index]))
        return
    end
    qs = qubits[active]
    rotations = Op[]
    for index in active
        terms[index] == "X" && push!(rotations, Op("H", qubits[index]))
        terms[index] == "Y" && push!(rotations, Op("HSP", qubits[index]))
    end
    append!(all_ops, rotations)
    ladder = [Op("CX", qs[i], qs[i+1]) for i in 1:length(qs)-1]
    append!(all_ops, ladder)
    push!(all_ops, Op("RZ($angle)", _RZ(angle), last(qs)))
    append!(all_ops, reverse(ladder))
    append!(all_ops, adjoint(rotations))
end

function _get_qubit_index(rows::Int, cols::Int, row::Int, col::Int, style::Symbol)
    style == :rowwise && return (row-1)*cols + col
    style == :colwise && return (col-1)*rows + row
    throw(ArgumentError("enumeration_style must be :rowwise or :colwise"))
end

_create_lattice_map(rows::Int, cols::Int, style::Symbol=:rowwise) =
    [[_get_qubit_index(rows,cols,r,c,style) for c in 1:cols] for r in 1:rows]

# Pauli strings commute iff the number of anticommuting shared sites is even.
function _trotter_commutes(a, b)
    odd = false
    for (pa, qa) in zip(a[2], a[3]), (pb, qb) in zip(b[2], b[3])
        if qa == qb && pa != "I" && pb != "I" && pa != pb
            odd = !odd
        end
    end
    return !odd
end

function _apply_trotter_group!(result, local_terms, group, duration)
    for index in group
        coupling, strings, qs = local_terms[index]
        _apply_term(result, strings, 2 * duration * coupling, qs)
    end
end

"""
    hamiltonian_exp(N, total_time, terms; dt=0.01, order=1, full=true)
    hamiltonian_exp((rows, cols), total_time, terms; dt=0.01, order=1, full=true, enumeration_style=:rowwise)

Construct a product-formula approximation to `exp(-im*total_time*H)` for
`H = hamiltonian(N, terms)` with open boundaries. `terms` alternates real
couplings and comma-separated Pauli strings, for example `[0.7, "X,Z", 0.2, "Y"]`.
Identity factors preserve the separation between the active sites. In 2D,
strings run along both lattice directions; single-site terms are counted once.

Choose `order=1` (Lie–Trotter, the default), `order=2` (symmetric Strang), or
`order=4` (Yoshida's three-stage composition). At fixed total time their errors
scale as `dt`, `dt^2`, and `dt^4` in the small-step limit. Fourth order uses
`S2(w*dt) S2((1-2w)*dt) S2(w*dt)`, where `w = 1/(2-cbrt(2))`; its middle
stage runs backwards. This uses fewer stages than the five-stage Suzuki formula,
but is not guaranteed to minimize gates at a prescribed error for every model.
See [Yoshida (1990)](https://doi.org/10.1016/0375-9601(90)90092-3).

Consecutive commuting local terms form exact blocks, and adjacent repetitions
of a block are merged across stages and time steps. Zero couplings are omitted;
single-site factors use one RX/RY/RZ gate. General Pauli strings use basis
rotations and `2(k-1)` CNOTs for `k` nonidentity sites. Identity factors retain
their global phase. The returned operation list does not store a register size;
apply it to an explicitly sized state when some qubits are idle.

The last step is shortened to reach the requested time exactly. `full=false`
returns just one step of duration `dt` (with the sign of `total_time`, positive
when it is zero). For example, `hamiltonian_exp(4, 0.6, [-1.0, "Z,Z", -0.8, "X"];
dt=0.05, order=2)` constructs second-order open-chain Ising evolution.
"""
function hamiltonian_exp(N::Int, total_time::Real, terms::Vector; dt=0.01, order=1, full=true)
    return hamiltonian_exp((1,N), total_time, terms; dt=dt, order=order, full=full)
end

function hamiltonian_exp(shape::Union{Tuple{Int,Int},Vector{Int}}, total_time::Real, terms::Vector;
                         dt=0.01, order=1, full=true, enumeration_style::Symbol=:rowwise)
    length(shape) == 2 || throw(ArgumentError("lattice dimensions must contain two integers"))
    rows, cols = shape
    rows > 0 && cols > 0 || throw(ArgumentError("lattice dimensions must be positive"))
    dt isa Real && isfinite(dt) && dt > 0 && isfinite(total_time) || throw(ArgumentError("time must be finite and dt must be positive"))
    order isa Integer && order in (1, 2, 4) || throw(ArgumentError("order must be 1, 2, or 4"))
    iseven(length(terms)) || throw(ArgumentError("terms must alternate couplings and Pauli strings"))
    _get_qubit_index(rows,cols,1,1,enumeration_style)
    couplings = terms[1:2:end]
    all(x -> x isa Real && isfinite(x), couplings) || throw(ArgumentError("couplings must be finite real numbers"))
    all(x -> x isa AbstractString, terms[2:2:end]) || throw(ArgumentError("Pauli terms must be strings"))
    paulis = [String.(strip.(split(uppercase(term), ','))) for term in terms[2:2:end]]
    all(p -> all(in(("I", "X", "Y", "Z")), p), paulis) || throw(ArgumentError("Trotter terms must contain I, X, Y, or Z"))

    # Place terms once, preserving the same order and boundaries as hamiltonian.
    local_terms = Tuple{Real,Vector{String},Vector{Int}}[]
    for (coupling, strings) in zip(couplings,paulis)
        iszero(coupling) && continue
        span = length(strings)
        for row in 1:rows, col in 1:cols-span+1
            qs = [_get_qubit_index(rows,cols,row,c,enumeration_style) for c in col:col+span-1]
            push!(local_terms, (coupling, strings, qs))
        end
        if span > 1
            for col in 1:cols, row in 1:rows-span+1
                qs = [_get_qubit_index(rows,cols,r,col,enumeration_style) for r in row:row+span-1]
                push!(local_terms, (coupling, strings, qs))
            end
        end
    end

    result = QuantumOps[]
    (isempty(local_terms) || (full && iszero(total_time))) && return result
    # A translated mixed string can fail to commute with its neighbours.
    # Check every member before putting it in an exactly exponentiable block.
    groups = [Int[]]
    for index in eachindex(local_terms)
        if !all(j -> _trotter_commutes(local_terms[index], local_terms[j]), last(groups))
            push!(groups, Int[])
        end
        push!(last(groups), index)
    end

    direction = total_time < 0 ? -1 : 1
    duration = full ? total_time : direction * dt
    if length(groups) == 1
        _apply_trotter_group!(result, local_terms, only(groups), duration)
        return result
    end

    # A symmetric sweep has a full central block and half outer blocks.
    n = length(groups)
    sweep = order == 1 ? [(i, 1.0) for i in 1:n] :
        vcat([(i, 0.5) for i in 1:n-1], [(n, 1.0)], [(i, 0.5) for i in n-1:-1:1])
    w = 1 / (2 - cbrt(2))
    stages = order == 4 ? (w, 1-2w, w) : (1.0,)
    count = full ? floor(Int, abs(total_time)/dt) : 1
    remainder = full ? max(zero(duration), abs(total_time) - count*dt) : zero(duration)
    # Retain every positive remainder, including times smaller than eps(dt).
    pending_group = 0
    pending_duration = zero(float(duration))
    for step in 1:count+!iszero(remainder)
        width = direction * (step <= count ? dt : remainder)
        for stage in stages, (group, weight) in sweep
            piece = width * stage * weight
            if group == pending_group
                pending_duration += piece
            else
                pending_group == 0 || _apply_trotter_group!(result, local_terms, groups[pending_group], pending_duration)
                pending_group = group
                pending_duration = piece
            end
        end
    end
    pending_group == 0 || _apply_trotter_group!(result, local_terms, groups[pending_group], pending_duration)
    return result
end

# Historical experimental name now shares the verified implementation.
hamiltonian_exp_test(args...; kwargs...) = hamiltonian_exp(args...; kwargs...)
