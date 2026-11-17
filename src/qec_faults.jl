# Pauli faults are represented by binary (x,z) vectors. Global phases cannot
# affect a syndrome or the logical fault and are deliberately discarded here.
_qec_fault_ops(errors::Op) = (errors,)
_qec_fault_ops(errors::AbstractVector) = errors

function _qec_pauli_vector(errors, n::Integer)
    n >= 0 || throw(ArgumentError("nqubits must be nonnegative"))
    v = zeros(Int, 2n)
    for op in _qec_fault_ops(errors)
        op isa Op && op.q == 1 && op.control == -2 && op.target_qubit == -1 &&
            op.name in ("I", "X", "Y", "Z") && op.mat isa AbstractMatrix &&
            isapprox(op.mat, gates(op.name); atol=1e-12, rtol=1e-12) ||
            throw(ArgumentError("faults must be single-qubit I, X, Y, or Z operations without extra controls"))
        1 <= op.qubit <= n || throw(ArgumentError("fault qubits must be between 1 and $n"))
        op.name in ("X", "Y") && (v[op.qubit] ⊻= 1)
        op.name in ("Z", "Y") && (v[n+op.qubit] ⊻= 1)
    end
    return v
end

function _qec_pauli_syndrome(code::StabilizerCode, v::AbstractVector)
    n = code.n
    return mod.(code.generator[:, 1:n] * v[n+1:end] +
                code.generator[:, n+1:end] * v[1:n], 2)
end

"""
    pauli_syndrome(code::StabilizerCode, errors) -> Vector{Int}

Compute the syndrome of a physical Pauli fault without simulating a state.
`errors` is one single-qubit Pauli `Op` or a vector of them; repeated operations
on a qubit are multiplied, ignoring global phase. The bits follow the input
stabilizer order: zero means commute and one means anticommute. A zero syndrome
does not establish that the fault acts trivially on the logical state.
"""
function pauli_syndrome(code::StabilizerCode, errors::Union{Op,AbstractVector})
    return _qec_pauli_syndrome(code, _qec_pauli_vector(errors, code.n))
end

function _qec_check_clifford(op)
    op isa Op || throw(ArgumentError("Pauli propagation requires unitary Clifford operations"))
    single = op.name in ("I", "X", "Y", "Z", "H", "S", "SD")
    pair = op.name in ("CX", "CNOT", "CY", "CZ", "SWAP", "ISWAP")
    (single || pair) && op.q == (single ? 1 : 2) && op.control == -2 &&
        op.mat isa AbstractMatrix &&
        isapprox(op.mat, gates(op.name); atol=1e-12, rtol=1e-12) ||
        throw(ArgumentError("unsupported Clifford operation: $(op.name)"))
    op.qubit >= 1 && (single ? op.target_qubit == -1 :
        op.target_qubit >= 1 && op.target_qubit != op.qubit) ||
        throw(ArgumentError("Clifford operation has invalid qubit indices"))
    return nothing
end

"""
    propagate_pauli(errors, circuit; nqubits=nothing) -> Vector{Op}

Propagate a Pauli fault through a Clifford circuit: if the gates implement `U`,
the output represents `U * E * U'`, up to an irrelevant global phase. Operations
are applied in their listed order. To inject a fault after gate `j`, propagate
through `circuit[j+1:end]`.

Supported gates are `I`, `X`, `Y`, `Z`, `H`, `S`, `SD`, `CX`/`CNOT`, `CY`, `CZ`,
`SWAP`, and `ISWAP`, with their built-in matrices and no extra controls.
Measurements, resets, channels, custom matrices, and non-Clifford gates (such as
√iSWAP) are rejected. The `noisy` flags do not sample noise: this is exact
algebra for a specified Pauli fault and an otherwise ideal Clifford circuit.

The returned single-qubit Pauli operations have `noisy=false`, are sorted by
qubit, and omit identity factors. The register size is inferred from both inputs
unless `nqubits` is supplied; it must include every operated qubit.
"""
function propagate_pauli(errors::Union{Op,AbstractVector}, circuit::AbstractVector;
                         nqubits::Union{Nothing,Integer}=nothing)
    largest_qubit = 0
    for op in _qec_fault_ops(errors)
        op isa Op || throw(ArgumentError("faults must be Pauli operations"))
        largest_qubit = max(largest_qubit, op.qubit)
    end
    for op in circuit
        _qec_check_clifford(op)
        largest_qubit = max(largest_qubit, op.qubit, op.target_qubit)
    end
    n = isnothing(nqubits) ? largest_qubit : nqubits
    n >= largest_qubit || throw(ArgumentError("nqubits must include every operated qubit"))
    v = _qec_pauli_vector(errors, n)
    for op in circuit
        q, t = op.qubit, op.target_qubit
        if op.name == "SWAP"
            v[q], v[t] = v[t], v[q]
            v[n+q], v[n+t] = v[n+t], v[n+q]
        else
            _qec_conjugate_gate!(v, op.name, q, t)
        end
    end
    return Op[Op(name, q; noisy=false) for (name, q) in logical_vec_2_ops(v)]
end

"""
    logical_pauli(code::StabilizerCode, residual) -> (x=bits, z=bits)

Classify a zero-syndrome physical Pauli in the encoder's logical basis. For
logical qubit `j`, `(x[j], z[j])` denotes `I=(0,0)`, `X=(1,0)`, `Z=(0,1)`, or
`Y=(1,1)`, ignoring global phase and stabilizers. Nonzero syndrome raises an
`ArgumentError`: first combine the fault with its candidate Pauli recovery.

This uses the logical representatives derived from the encoder, even when
`code.logicals` contains user-supplied gate mappings. A zero result means the
residual acts trivially within this code space. Validating a noisy correction
protocol also requires checking faults in its preparation, gates, and readout.
"""
function logical_pauli(code::StabilizerCode, residual::Union{Op,AbstractVector})
    v = _qec_pauli_vector(residual, code.n)
    all(iszero, _qec_pauli_syndrome(code, v)) ||
        throw(ArgumentError("logical Pauli classification requires zero syndrome"))
    n = code.n
    X, Z = code.info.logical_XZ_vecs
    x = mod.(Z[:, 1:n] * v[n+1:end] + Z[:, n+1:end] * v[1:n], 2)
    z = mod.(X[:, 1:n] * v[n+1:end] + X[:, n+1:end] * v[1:n], 2)
    return (x=Vector{Int}(x), z=Vector{Int}(z))
end
