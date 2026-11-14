"""
    qec_round(code::StabilizerCode, state; decoder=nothing)

Measure the code's stabilizers ideally, apply a Pauli recovery, and return
`(state=corrected_data, syndrome=bits)`. The returned state contains exactly
`code.n` data qubits, so rounds can be interleaved with logical gates without
accumulating ancillas. Syndrome bits follow the input-generator order, with
zero for eigenvalue +1 and one for eigenvalue -1.

Statevectors and MPS sample a selective syndrome and return a normalized pure
state. MPS projection stays in tensor form without intentional truncation.
Density matrices sum the recovered syndrome sectors deterministically and
return `syndrome=nothing`, since an ensemble has no single measurement record.

The default decoder searches for one physical Pauli error. A custom
`decoder(bits)` returns a Pauli `Op` or vector of them, as for `code.correct`.
Syndrome consistency alone does not rule out a residual logical error.

This is an **ideal recovery operation**, implemented with projectors
`(I + (-1)^bit S)/2`. It has no gate, measurement, or recovery noise and is not
a fault-tolerant syndrome-extraction circuit. Use it as an ideal reference
between noisy gate gadgets; it cannot establish an end-to-end threshold.
"""
function qec_round(code::StabilizerCode, state::Union{AbstractVectorS,it.MPS}; decoder=nothing)
    get_N(state) == code.n || throw(DimensionMismatch("an ideal QEC round requires $(code.n) data qubits"))
    state_norm = la.norm(state)
    isfinite(state_norm) && state_norm > 0 || throw(ArgumentError("the input state must have a finite, positive norm"))
    psi = la.normalize(state)
    bits = Int[]
    for row in eachrow(code.generator)
        # S is Hermitian, including Y factors; no binary-Pauli phase is omitted.
        pauli_ops = Op[Op(name, q; noisy=false) for (name, q) in logical_vec_2_ops(collect(row))]
        Spsi = apply(psi, pauli_ops)
        p_zero = clamp((1 + real(inner(psi, Spsi))) / 2, 0, 1)
        bit = Int(rand() >= p_zero)
        psi = _qec_project_pure(psi, Spsi, 1 - 2bit)
        push!(bits, bit)
    end
    return (state=code.correct(psi, bits; decoder=decoder), syndrome=bits)
end

_qec_project_pure(psi::AbstractVector, Spsi::AbstractVector, sign) = la.normalize(psi + sign * Spsi)

function _qec_project_pure(psi::it.MPS, Spsi::it.MPS, sign)
    projected = +(psi, sign * Spsi; alg="directsum")
    # Compress the direct-sum representation, retaining all singular values.
    # No full statevector or exponentially sized Pauli matrix is constructed.
    it.truncate!(projected; cutoff=0.0)
    return la.normalize(projected)
end

function qec_round(code::StabilizerCode, rho::AbstractMatrixS; decoder=nothing)
    size(rho) == (2^code.n, 2^code.n) ||
        throw(DimensionMismatch("an ideal QEC round requires a $(code.n)-qubit density matrix"))
    trace_rho = real(la.tr(rho))
    isfinite(trace_rho) && trace_rho > 0 ||
        throw(ArgumentError("the input density matrix must have a finite, positive trace"))

    # Keep the density calculation in the data Hilbert space. Appending all
    # syndrome ancillas would multiply its dimension by 2^code.m needlessly.
    sectors = Tuple{AbstractMatrix,Vector{Int}}[(sa.sparse(rho), Int[])]
    identity = sa.sparse(1.0 * la.I, 2^code.n, 2^code.n)
    for stabilizer in code.stabilizers
        S = string_to_matrix(stabilizer)
        next_sectors = typeof(sectors)()
        for (branch, bits) in sectors, bit in 0:1
            projector = (identity + (1-2bit) * S) / 2
            projected = projector * branch * projector'
            # Discard floating-point remnants of zero-probability sectors.
            real(la.tr(projected)) > 100eps(Float64) * trace_rho || continue
            push!(next_sectors, (projected, [bits; bit]))
        end
        sectors = next_sectors
    end
    recovered = zero(sa.sparse(rho))
    for (branch, bits) in sectors
        recovered += code.correct(branch, bits; decoder=decoder)
    end
    return (state=recovered, syndrome=nothing)
end
