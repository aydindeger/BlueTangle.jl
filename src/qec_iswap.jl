# A fault-tolerant logical iSWAP: two [[6,2,2]] blocks (Knill's C6 code) and transversal iSWAP gates.
#
# One block occupies a row of six qubits. Along the row the positions hold a, b, a+b, a, a+b, b in the
# codewords of logical |00>: sum over bits a, b of |a, b, a+b, a, a+b, b>. The X and Z checks have the same
# supports, so the code is self-dual. Logical qubit i of a block has X_L = X and Z_L = Z on the odd-weight
# support _C6_LOGICAL[i]. Odd weight is what makes transversal S act as logical S† and transversal iSWAP between
# two blocks act as logical iSWAP† on corresponding logical qubits.

const _C6_CHECKS = ([1, 0, 1, 1, 1, 0], [0, 1, 1, 0, 1, 1])
const _C6_LOGICAL = ([1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1])
const _C6_INFO = (1, 3)
# Nearest-neighbour CNOT layers (control, target) preparing logical |00> from |0...0> after H on _C6_INFO.
# No single fault in them leaves an undetectable logical error; seven CNOTs suffice without that property.
const _C6_ENCODER = (((1, 2), (3, 4)), ((2, 3), (4, 5)), ((3, 4), (5, 6)), ((4, 5), (3, 2)))

_c6_wires(block::Integer) = 6(block - 1) .+ (1:6)

function _c6_check_block(block)
    block isa Integer && block >= 1 || throw(ArgumentError("blocks are numbered from one"))
end

function _c6_check_basis(basis)
    basis in (:Z, :X, :Y) || throw(ArgumentError("basis must be :Z, :X or :Y"))
end

"""
    c6_code(blocks::Int=1)

Return the stabilizer group of `blocks` independent [[6,2,2]] blocks as a `StabilizerCode`. Block `b` occupies
physical wires `6(b-1)+1:6b`. Each block has X and Z checks on the same two supports, `X,I,X,X,X,I` and
`I,X,X,I,X,X`, in that order: X then Z for the first support, then for the second.

The code detects every single-qubit error (`d=2`) and encodes two logical qubits per block. The fault-tolerant
gadgets `c6_prepare`, `c6_iswap`, `c6_readout` and `c6_decode` use the logical basis of `c6_logicals`: for
logical qubit `i` of a block, X and Z on positions 1–3 (`i=1`) or 4–6 (`i=2`). The returned object's own encoder
(`code.encode`) synthesizes a different logical basis. Its syndromes (`pauli_syndrome`), `code_distance`, and
the test of whether `logical_pauli` is zero do not depend on that choice.
"""
function c6_code(blocks::Int=1)
    blocks > 0 || throw(ArgumentError("the number of blocks must be positive"))
    rows = String[]
    for b in 1:blocks, check in _C6_CHECKS, name in ("X", "Z")
        row = fill("I", 6blocks)
        row[_c6_wires(b)[check .== 1]] .= name
        push!(rows, join(row, ","))
    end
    return StabilizerCode(rows; d=2)
end

"""
    c6_logicals(blocks::Int=1) -> Dict

Logical Pauli representatives used by the C6 gadgets: `("X", j)` and `("Z", j)` for logical qubit
`j = 2(b-1)+i`, as vectors of physical `Op`s on block `b`. Logical Y is `-Y⊗Y⊗Y` on the same support.
"""
function c6_logicals(blocks::Int=1)
    blocks > 0 || throw(ArgumentError("the number of blocks must be positive"))
    out = Dict{Tuple{String,Int},Vector{Op}}()
    for b in 1:blocks, i in 1:2, name in ("X", "Z")
        out[(name, 2(b-1)+i)] = [Op(name, q; noisy=false) for q in _c6_wires(b)[_C6_LOGICAL[i] .== 1]]
    end
    return out
end

"""
    c6_prepare(basis::Symbol=:Z; bits=(0, 0), block::Int=1) -> Vector{Op}

Prepare both logical qubits of a [[6,2,2]] block in an eigenstate of logical `basis` (`:Z`, `:X` or `:Y`), with
eigenvalues `(-1)^bits[1]` and `(-1)^bits[2]`. Start from |0⟩ on the block's six wires.

The circuit is fault tolerant without a flag qubit. H on positions 1 and 3 and eight nearest-neighbour CNOTs in
four layers prepare logical |00⟩. Transversal H gives logical H on both qubits. Transversal SD then gives logical
S on both. X (for `:Z`) or Z (for `:X`, `:Y`) on a logical support sets the value. Any single fault in these
gates either flips a parity check of the final transversal readout or leaves the logical state unchanged. The
Z flips are marked `noisy=false`: superconducting processors apply Z rotations in software.
"""
function c6_prepare(basis::Symbol=:Z; bits=(0, 0), block::Int=1)
    _c6_check_basis(basis)
    _c6_check_block(block)
    length(bits) == 2 && all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be two zeros or ones"))
    w = _c6_wires(block)
    ops = Op[Op("H", w[q]) for q in _C6_INFO]
    for layer in _C6_ENCODER, (c, t) in layer
        push!(ops, Op("CX", w[c], w[t]))
    end
    basis in (:X, :Y) && append!(ops, [Op("H", q) for q in w])
    basis == :Y && append!(ops, [Op("SD", q) for q in w])
    for i in 1:2
        bits[i] == 1 || continue
        support = w[_C6_LOGICAL[i] .== 1]
        append!(ops, basis == :Z ? [Op("X", q) for q in support] : [Op("Z", q; noisy=false) for q in support])
    end
    return ops
end

"""
    c6_iswap(; blocks=(1, 2), dagger=false) -> Vector{Op}

Logical iSWAP between corresponding logical qubits of two [[6,2,2]] blocks: on (1st of `blocks[1]`, 1st of
`blocks[2]`) and, simultaneously, on their second logical qubits.

The physical gates are six `ISWAP`s, one on each pair of corresponding wires. On the code space that transversal
gate is logical iSWAP† on both pairs. Z on every wire of the first block, before and after, turns it into
logical iSWAP. With `dagger=true` the Z layers are omitted, giving logical iSWAP†. Each physical gate acts on one
wire of each block. A single fault therefore leaves at most one error per block, which the code detects.

Transversal √iSWAP is not a logical gate of this code: it leaves the code space. A physical iSWAP is two native
√iSWAPs on superconducting processors. The Z gates are marked `noisy=false` (applied in software).
"""
function c6_iswap(; blocks::Tuple{Int,Int}=(1, 2), dagger::Bool=false)
    _c6_check_block(blocks[1])
    _c6_check_block(blocks[2])
    blocks[1] != blocks[2] || throw(ArgumentError("the logical iSWAP acts between two different blocks"))
    a, b = _c6_wires(blocks[1]), _c6_wires(blocks[2])
    pairs = Op[Op("ISWAP", a[j], b[j]) for j in 1:6]
    dagger && return pairs
    z = Op[Op("Z", q; noisy=false) for q in a]
    return vcat(z, pairs, z)
end

"""
    c6_readout(basis::Symbol; block::Int=1) -> Vector{Op}

Single-qubit rotations that turn a transversal measurement of `basis` on the block into a Z-basis measurement:
none for `:Z`, H for `:X`, SD then H for `:Y`. Measure every wire of the block afterwards and pass the bits to
`c6_decode`.
"""
function c6_readout(basis::Symbol; block::Int=1)
    _c6_check_basis(basis)
    _c6_check_block(block)
    w = _c6_wires(block)
    basis == :Z && return Op[]
    basis == :X && return Op[Op("H", q) for q in w]
    return vcat(Op[Op("SD", q) for q in w], Op[Op("H", q) for q in w])
end

"""
    c6_decode(bits::AbstractVector{<:Integer}, bases) -> (accepted, logical)

Decode one shot of a transversal readout of `length(bases)` [[6,2,2]] blocks. `bits` holds the Z-basis outcomes
of all their wires, block 1 first, after `c6_readout(bases[b]; block=b)`. A shot is accepted when the two parity
checks of every block are even. `logical` holds two bits per block: the parity over each logical support. For Y
one is added, since logical Y is `-Y⊗Y⊗Y`. A single flipped bit is always rejected.
"""
function c6_decode(bits::AbstractVector{<:Integer}, bases)
    length(bits) == 6length(bases) || throw(DimensionMismatch("expected six bits per block"))
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be zeros or ones"))
    accepted = true
    logical = Int[]
    for (b, basis) in enumerate(bases)
        _c6_check_basis(basis)
        block = bits[_c6_wires(b)]
        accepted &= all(iszero(sum(block .* check) % 2) for check in _C6_CHECKS)
        append!(logical, [(sum(block .* support) + (basis == :Y)) % 2 for support in _C6_LOGICAL])
    end
    return (accepted=accepted, logical=logical)
end

"""
    c6_single_faults(circuit; blocks::Int=2) -> NamedTuple

Insert every single fault into a C6 circuit, one at a time. The circuit must be Clifford, act on `6blocks` wires,
and end with `c6_readout` on every block, so that all wires are then measured in Z. A fault is a non-identity
Pauli after a gate on that gate's wires (3 or 15 choices), an X on a fresh wire, or a flipped measurement.

Each fault is propagated to the measurement and classified:
- `detected`: some parity check of some block flips, so the shot is discarded;
- `harmless`: no check and no logical bit flips;
- `harmful`: no check flips but a logical bit does, an undetectable logical error.

Returns the counts and `harmful_faults`, the list of `(location, paulis)` for the harmful ones. Location 0 means
a fresh wire and `:measure` a flipped outcome. A fault-tolerant circuit has no harmful single fault.
"""
function c6_single_faults(circuit::AbstractVector; blocks::Int=2)
    n = 6blocks
    function classify(propagated)
        flips = zeros(Int, n)
        for op in propagated
            op.name in ("X", "Y") && (flips[op.qubit] ⊻= 1)
        end
        detected = wrong = false
        for b in 1:blocks
            block = flips[_c6_wires(b)]
            detected |= any(isodd(sum(block .* c)) for c in _C6_CHECKS)
            wrong |= any(isodd(sum(block .* s)) for s in _C6_LOGICAL)
        end
        return detected, wrong
    end
    counts = Dict(:detected => 0, :harmless => 0, :harmful => 0)
    harmful = Tuple{Any,Vector{Op}}[]
    record(location, fault, result) = begin
        detected, wrong = result
        key = detected ? :detected : wrong ? :harmful : :harmless
        counts[key] += 1
        key == :harmful && push!(harmful, (location, fault))
    end
    for q in 1:n
        fault = [Op("X", q; noisy=false)]
        record(0, fault, classify(propagate_pauli(fault, circuit; nqubits=n)))
        record(:measure, fault, classify(fault))
    end
    for (location, op) in enumerate(circuit)
        wires = op.q == 1 ? (op.qubit,) : (op.qubit, op.target_qubit)
        rest = circuit[location+1:end]
        for names in Iterators.product(ntuple(_ -> ("I", "X", "Y", "Z"), length(wires))...)
            all(==("I"), names) && continue
            fault = [Op(name, q; noisy=false) for (name, q) in zip(names, wires) if name != "I"]
            record(location, fault, classify(propagate_pauli(fault, rest; nqubits=n)))
        end
    end
    total = sum(values(counts))
    return (faults=total, detected=counts[:detected], harmless=counts[:harmless], harmful=counts[:harmful],
            harmful_faults=harmful)
end

# ----------------------------------------------------------------- Pauli frames
function _frame_gate!(x, z, name, q, t)
    if name == "H"
        x[q], z[q] = z[q], x[q]
    elseif name in ("S", "SD")
        z[q] .⊻= x[q]
    elseif name in ("CX", "CNOT")
        x[t] .⊻= x[q]
        z[q] .⊻= z[t]
    elseif name == "CZ"
        z[q] .⊻= x[t]
        z[t] .⊻= x[q]
    elseif name == "CY"
        z[t] .⊻= x[t]
        x[t] .⊻= x[q]
        z[q] .⊻= z[t]
        z[t] .⊻= x[t]
    elseif name == "SWAP"
        x[q], x[t] = x[t], x[q]
        z[q], z[t] = z[t], z[q]
    elseif name == "ISWAP"                  # iSWAP = SWAP · CZ · (S ⊗ S)
        z[q] .⊻= x[q]
        z[t] .⊻= x[t]
        z[q] .⊻= x[t]
        z[t] .⊻= x[q]
        x[q], x[t] = x[t], x[q]
        z[q], z[t] = z[t], z[q]
    end
    return nothing
end

"""
    sample_pauli_frames(circuit, p; nqubits, shots=1000, p2=p, prep=0.0, readout=0.0) -> (x, z)

Sample independent Pauli faults through a Clifford circuit and return the accumulated Pauli frame of every shot,
as `BitMatrix`es of size `shots × nqubits`. Bit `x[s, q]` means shot `s` ends with an X component on qubit `q`.
After a Z-basis measurement, `x` is the pattern of flipped outcomes.

Each operation with `noisy=true` is followed, with probability `p` (one-qubit gates) or `p2` (two-qubit gates),
by a uniformly random non-identity Pauli on its qubits: the depolarizing channel. Every qubit starts with an X
fault with probability `prep`, and its final outcome flips with probability `readout`. Global phases and Pauli
signs are discarded. The supported gates are those of `propagate_pauli`.

For Clifford circuits and Pauli noise this is exact and much faster than a state-vector simulation. Combine the
frame with the ideal outcome of a noiseless simulation, or test it directly, as `pauli_syndrome` and
`logical_pauli` do.
"""
function sample_pauli_frames(circuit::AbstractVector, p::Real; nqubits::Integer, shots::Integer=1000,
                             p2::Real=p, prep::Real=0.0, readout::Real=0.0)
    all(r -> isfinite(r) && 0 <= r <= 1, (p, p2, prep, readout)) ||
        throw(ArgumentError("error probabilities must lie in [0, 1]"))
    shots > 0 || throw(ArgumentError("shots must be positive"))
    for op in circuit
        _qec_check_clifford(op)
        max(op.qubit, op.target_qubit) <= nqubits || throw(ArgumentError("nqubits must include every operated qubit"))
    end
    x = [rand(shots) .< prep for _ in 1:nqubits]
    z = [falses(shots) for _ in 1:nqubits]
    for op in circuit
        _frame_gate!(x, z, op.name, op.qubit, op.target_qubit)
        op.noisy || continue
        qubits = op.q == 1 ? (op.qubit,) : (op.qubit, op.target_qubit)
        prob = op.q == 1 ? p : p2
        iszero(prob) && continue
        for s in findall(rand(shots) .< prob)
            choice = rand(1:4^length(qubits)-1)          # a uniformly random non-identity Pauli
            for (m, q) in enumerate(qubits)
                digit = (choice >> (2(m-1))) & 3          # 1 = X, 2 = Z, 3 = Y
                x[q][s] ⊻= isodd(digit)
                z[q][s] ⊻= digit >= 2
            end
        end
    end
    if readout > 0
        for q in 1:nqubits
            x[q] .⊻= rand(shots) .< readout
        end
    end
    return (x=reduce(hcat, x), z=reduce(hcat, z))
end
