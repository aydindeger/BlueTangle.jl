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
    trotter_ising(N, total_time, J, h; dt=0.1)

First-order evolution for `H = -J*sum(Z_i*Z_{i+1}) - h*sum(X_i)`.
"""
trotter_ising(N::Int, total_time::Real, J::Real, h::Real; dt=0.1) =
    hamiltonian_exp(N, total_time, [-h, "X", -J, "Z,Z"]; dt=dt)

# exp(-i*angle/2*P): rotate each nonidentity Pauli to Z, collect its parity,
# rotate the parity qubit, and undo the parity circuit and basis rotations.
function _apply_term(all_ops::Vector{QuantumOps}, terms::AbstractVector, angle::Real, qubits::Vector{Int})
    all(in(("I", "X", "Y", "Z")), terms) || throw(ArgumentError("Trotter terms must contain I, X, Y, or Z"))
    active = findall(!=("I"), terms)
    if isempty(active)
        push!(all_ops, Op("global phase", exp(-im*angle/2) * gate.I, first(qubits)))
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

"""
    hamiltonian_exp(N, total_time, terms; dt=0.01, full=true)
    hamiltonian_exp((rows, cols), total_time, terms; dt=0.01, full=true, enumeration_style=:rowwise)

Construct a first-order Trotter approximation to `exp(-im*total_time*H)` for
`H = hamiltonian(N, terms)` with open boundaries. `terms` alternates real
couplings and comma-separated Pauli strings, for example `[0.7, "X,Z", 0.2, "Y"]`.
Identity factors preserve the separation between the active sites. In 2D,
strings run along both lattice directions; single-site terms are counted once.

The last step is shortened to reach the requested time exactly. `full=false`
returns just one step of duration `dt` (with the sign of `total_time`).
"""
function hamiltonian_exp(N::Int, total_time::Real, terms::Vector; dt=0.01, full=true)
    return hamiltonian_exp((1,N), total_time, terms; dt=dt, full=full)
end

function hamiltonian_exp(shape::Union{Tuple{Int,Int},Vector{Int}}, total_time::Real, terms::Vector;
                         dt=0.01, full=true, enumeration_style::Symbol=:rowwise)
    rows, cols = shape
    rows > 0 && cols > 0 || throw(ArgumentError("lattice dimensions must be positive"))
    isfinite(dt) && dt > 0 && isfinite(total_time) || throw(ArgumentError("time must be finite and dt must be positive"))
    iseven(length(terms)) || throw(ArgumentError("terms must alternate couplings and Pauli strings"))
    _get_qubit_index(rows,cols,1,1,enumeration_style)
    couplings = terms[1:2:end]
    all(x -> x isa Real && isfinite(x), couplings) || throw(ArgumentError("couplings must be finite real numbers"))
    paulis = [strip.(split(uppercase(term), ',')) for term in terms[2:2:end]]
    steps = if full
        count = floor(Int, abs(total_time)/dt)
        widths = fill(Float64(dt), count)
        remainder = abs(total_time) - count*dt
        remainder > eps(max(abs(float(total_time)), float(dt))) && push!(widths, remainder)
        widths
    else
        [Float64(dt)]
    end
    direction = total_time < 0 ? -1 : 1
    result = QuantumOps[]
    for width in steps, (coupling, strings) in zip(couplings,paulis)
        span = length(strings)
        angle = 2 * direction * width * coupling
        for row in 1:rows, col in 1:cols-span+1
            qs = [_get_qubit_index(rows,cols,row,c,enumeration_style) for c in col:col+span-1]
            _apply_term(result, strings, angle, qs)
        end
        if span > 1
            for col in 1:cols, row in 1:rows-span+1
                qs = [_get_qubit_index(rows,cols,r,col,enumeration_style) for r in row:row+span-1]
                _apply_term(result, strings, angle, qs)
            end
        end
    end
    return result
end

# Historical experimental name now shares the verified implementation.
hamiltonian_exp_test(args...; kwargs...) = hamiltonian_exp(args...; kwargs...)
