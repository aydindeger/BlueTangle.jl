# The bounded CSS support search follows Aydin Deger's native CCQC.jl
# src/distance.jl, adapted here with the author's permission. No CCQC Python
# gatefinder, automorphism, or diagonal-gate-finding code is imported or copied.
# BlueTangle additionally handles general stabilizer codes with mixed Paulis.

# Enumerate supports without materializing all binomial(n, weight) choices.
function _qec_distance_supports(visit, n, weight)
    support = collect(1:weight)
    while true
        visit(support) || return false
        position = weight
        while position > 0 && support[position] == n - weight + position
            position -= 1
        end
        position == 0 && return true
        support[position] += 1
        for j in position+1:weight
            support[j] = support[j-1] + 1
        end
    end
end

# A commuting Pauli is a nontrivial logical exactly when it anticommutes with
# at least one member of a complete logical X/Z basis. This excludes stabilizers
# even for degenerate codes, without enumerating their exponentially large group.
function _qec_distance_candidate(checks, logicals, x, z)
    n = length(x)
    for row in axes(checks, 1)
        parity = false
        for q in 1:n
            parity ⊻= (checks[row, q] == 1 && z[q]) ⊻ (checks[row, n+q] == 1 && x[q])
        end
        parity && return false
    end
    for row in axes(logicals, 1)
        parity = false
        for q in 1:n
            parity ⊻= (logicals[row, q] == 1 && z[q]) ⊻ (logicals[row, n+q] == 1 && x[q])
        end
        parity && return true
    end
    return false
end

"""
    code_distance(code::StabilizerCode; max_weight=4, max_candidates=250_000)

Search increasing Pauli weights for a nontrivial logical operator. A candidate
must commute with every stabilizer and anticommute with a stored logical X or
Z. The result therefore excludes low-weight stabilizers of degenerate codes.
The calculation is independent of the stored `code.d` value.

For CSS generators, pure X and pure Z supports suffice and their separate
distances `dx` and `dz` are reported. For other stabilizer codes, all X/Y/Z
assignments on each support are tested. The work is bounded by `max_weight`
and `max_candidates`; each tested Pauli counts as one candidate.

Return a named tuple with `distance`, `exact`, `lower_bound`, and `witness`
(physical Pauli operations). A found minimum is exact. If no witness is found,
`distance` and `witness` are `nothing`, and only the lower bound proved by fully
completed weights is returned. `searched_through` is the largest fully searched
weight; `truncated` means the candidate budget interrupted the search.
`stop_reason` is `:found`, `:max_weight`, or `:max_candidates`.

For CSS codes, `dx`, `dz`, `dx_lower_bound`, `dz_lower_bound`, `x_witness`, and
`z_witness` give the same information for each Pauli type; the two witnesses
are one-indexed supports. These fields are `nothing` for non-CSS codes.
The method is `:css_support` or `:pauli_support`. The number of candidates can
grow exponentially with code size.
"""
function code_distance(code::StabilizerCode; max_weight::Int=4, max_candidates::Int=250_000)
    code.k > 0 || throw(ArgumentError("distance is undefined for a code encoding no logical qubits"))
    max_weight >= 1 || throw(ArgumentError("max_weight must be positive"))
    max_candidates >= 1 || throw(ArgumentError("max_candidates must be positive"))
    n = code.n
    checks = Matrix{Int}(code.generator)
    logicals = Matrix{Int}(vcat(code.info.logical_XZ_vecs...))
    css = all(row -> all(iszero, row[1:n]) || all(iszero, row[n+1:end]), eachrow(checks))
    x, z = falses(n), falses(n)
    dx = dz = nothing
    x_witness = z_witness = witness = nothing
    distance = nothing
    candidates = searched_through = 0
    truncated = false

    for weight in 1:min(max_weight, n)
        completed = _qec_distance_supports(n, weight) do support
            assignments = css ? (("X",), ("Z",)) :
                Iterators.product(ntuple(_ -> ("X", "Y", "Z"), weight)...)
            for names in assignments
                css && ((names[1] == "X" && dx !== nothing) ||
                        (names[1] == "Z" && dz !== nothing)) && continue
                if candidates == max_candidates
                    truncated = true
                    return false
                end
                candidates += 1
                fill!(x, false); fill!(z, false)
                for (i, q) in enumerate(support)
                    name = css ? names[1] : names[i]
                    x[q] = name in ("X", "Y")
                    z[q] = name in ("Y", "Z")
                end
                if _qec_distance_candidate(checks, logicals, x, z)
                    if distance === nothing
                        distance = weight
                        witness = [Op(x[q] ? (z[q] ? "Y" : "X") : "Z", q; noisy=false) for q in support]
                    end
                    css || return false
                    if names[1] == "X"
                        dx, x_witness = weight, copy(support)
                    else
                        dz, z_witness = weight, copy(support)
                    end
                    dx !== nothing && dz !== nothing && return false
                end
            end
            return true
        end
        completed && (searched_through = weight)
        (truncated || (!css && distance !== nothing) || (dx !== nothing && dz !== nothing)) && break
    end
    frontier = searched_through + 1
    return (distance=distance, exact=distance !== nothing,
            lower_bound=distance === nothing ? frontier : distance, witness=witness,
            dx=dx, dz=dz,
            dx_lower_bound=css ? (dx === nothing ? frontier : dx) : nothing,
            dz_lower_bound=css ? (dz === nothing ? frontier : dz) : nothing,
            x_witness=x_witness, z_witness=z_witness,
            searched_through=searched_through, candidates=candidates, truncated=truncated,
            stop_reason=truncated ? :max_candidates :
                ((!css && distance !== nothing) || (dx !== nothing && dz !== nothing) ? :found : :max_weight),
            method=css ? :css_support : :pauli_support)
end
