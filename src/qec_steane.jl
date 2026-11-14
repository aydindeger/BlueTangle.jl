# The row order is part of the decoder convention: three X checks, then three
# Z checks, independently for each contiguous seven-qubit block.
const _STEANE_CHECKS = [
    "X,X,X,X,I,I,I", "X,X,I,I,X,X,I", "X,I,X,I,X,I,X",
    "Z,Z,Z,Z,I,I,I", "Z,Z,I,I,Z,Z,I", "Z,I,Z,I,Z,I,Z",
]

function _steane_checks(blocks::Int)
    return [join(vcat(fill("I", 7(b-1)), split(row, ","),
                      fill("I", 7(blocks-b))), ",")
            for b in 1:blocks for row in _STEANE_CHECKS]
end

"""
    steane_code(blocks::Int=1)

Construct `blocks` independent [[7,1,3]] Steane code blocks. Logical qubit `b`
occupies physical wires `7(b-1)+1:7b`; the input and output logical states use
this same block order. This is a [[7blocks, blocks, 3]] code, not concatenation.

The existing `code.ops` and `code.apply` interfaces support logical I, X, Y, Z,
H, S, SD, CX/CNOT, CZ and SWAP. These use single-qubit layers, or two-qubit gates
connecting corresponding positions in two different blocks. Logical S uses
physical SD on each wire (and conversely); Y uses the sequence SD, X, S.

Use `decoder = steane_decoder(code)` with `qec_round` for ideal recovery of
one arbitrary Pauli error per block. The default generic decoder instead
searches for a single error in the entire register.

Transversal gates contain a single gate fault to at most one qubit per block.
The synthesized encoder and `code.ops_syndrome` are NOT fault-tolerant gadgets;
`qec_round` assumes perfect syndrome measurement and recovery. No logical T gate
or magic-state preparation is supplied.
"""
function steane_code(blocks::Int=1)
    blocks > 0 || throw(ArgumentError("the number of Steane blocks must be positive"))
    logicals = Dict{Any,Any}()
    for b in 1:blocks
        wires = 7(b-1)+1:7b
        for name in ("I", "X", "Z", "H", "S", "SD")
            physical = name == "S" ? "SD" : name == "SD" ? "S" : name
            logicals[name, b] = name == "I" ? Op[] : [Op(physical, q) for q in wires]
        end
        # Y = i X Z. Seven physical Y's implement -Y in this convention, so
        # keep the exact logical phase by composing S†, X, S transversally.
        logicals["Y", b] = vcat(logicals["SD", b], logicals["X", b], logicals["S", b])
    end
    for b in 1:blocks, c in 1:blocks
        b == c && continue
        for name in ("CX", "CNOT", "CZ", "SWAP")
            logicals[name, b, c] = [Op(name, 7(b-1)+q, 7(c-1)+q) for q in 1:7]
        end
    end
    return StabilizerCode(_steane_checks(blocks), logicals; d=3)
end

"""
    steane_decoder(code::StabilizerCode)

Return a syndrome-to-Pauli callback for the check order used by `steane_code`.
Pass it as `decoder=steane_decoder(code)` to `qec_round`, `code.correct` or
`code.syndrome`. The three Z checks locate X errors and the three X checks
locate Z errors, independently in each block. A Y error has both components.

Every single-qubit Pauli per block is corrected. Higher-weight errors can have
the same syndrome and cause a logical error after this recovery. A zero
syndrome is not a certificate of a correct logical state.
"""
function steane_decoder(code::StabilizerCode)
    blocks = code.k
    code.n == 7blocks && code.m == 6blocks && blocks > 0 ||
        throw(ArgumentError("steane_decoder requires independent seven-qubit Steane blocks"))
    expected = stabilizers_to_generator(_steane_checks(blocks))
    code.generator == expected || throw(ArgumentError("the checks must have the order and wire layout of steane_code"))
    checks = expected[1:3, 1:7]

    function decode_syndrome(bits::AbstractVector)
        length(bits) == code.m || throw(DimensionMismatch("expected $(code.m) syndrome bits"))
        all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("syndrome bits must be zero or one"))
        recovery = Op[]
        for block in 1:blocks
            offset = 6(block-1)
            sx, sz = bits[offset+1:offset+3], bits[offset+4:offset+6]
            xsite = all(iszero, sz) ? 0 : findfirst(q -> checks[:, q] == sz, 1:7)
            zsite = all(iszero, sx) ? 0 : findfirst(q -> checks[:, q] == sx, 1:7)
            for q in 1:7
                x, z = q == xsite, q == zsite
                if x || z
                    push!(recovery, Op(x ? (z ? "Y" : "X") : "Z", 7(block-1)+q; noisy=false))
                end
            end
        end
        return recovery
    end
    return decode_syndrome
end
