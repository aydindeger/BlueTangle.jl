"""
    shor_code()

Construct Shor's [[9,1,3]] stabilizer code: three repetition blocks detect bit
flips, and two X checks compare the phases of adjacent triples. Logical gates
and the logical input basis are those derived by `StabilizerCode`'s encoder.

This nine-qubit code is distinct from the Shor *syndrome-extraction method*:
`shor_qec_round` currently implements that method for Steane blocks only.
The stored `d=3` is metadata; `code_distance(shor_code())` verifies it and
returns a physical logical-Pauli witness.
"""
function shor_code()
    return StabilizerCode([
        "Z,Z,I,I,I,I,I,I,I", "I,Z,Z,I,I,I,I,I,I",
        "I,I,I,Z,Z,I,I,I,I", "I,I,I,I,Z,Z,I,I,I",
        "I,I,I,I,I,I,Z,Z,I", "I,I,I,I,I,I,I,Z,Z",
        "X,X,X,X,X,X,I,I,I", "I,I,I,X,X,X,X,X,X",
    ]; d=3)
end

"""
    reed_muller_code()

Construct the [[15,1,3]] punctured quantum Reed–Muller CSS code used by the
15-to-1 magic-state protocol. Physical columns are the nonzero four-bit words
`1:15`, most-significant bit first. The four coordinate rows are X checks;
the ten Z checks span the binary kernel of these rows and the all-ones row.

Logical X and Z are equivalent to X and Z on all fifteen wires. Physical
`T⊗15` implements logical T†, so `code.ops("T", 1)` uses physical TD and
`code.ops("TD", 1)` uses physical T. This follows from the weights 0/8 in
logical zero and 15/7 in logical one. No diagonal-gate finder is used.

`code_distance(code; max_weight=7)` verifies `dx=7`, `dz=3`, hence `d=3`.
The general synthesized encoder is not certified fault tolerant.
"""
function reed_muller_code()
    coordinates = [Int((column >> (4-row)) & 1) for row in 1:4, column in 1:15]
    G = vcat(ones(Int, 1, 15), coordinates)
    standard, rank, _, permutation = reduced_row_echelon(copy(G))
    # In permuted coordinates, G=[I B], so columns of [B;I] form its GF(2)
    # kernel. Map them back to the stated physical column order.
    kernel = permutation * vcat(standard[:, rank+1:end], Matrix{Int}(la.I, 15-rank, 15-rank))
    zchecks = permutedims(mod.(kernel, 2))
    rows = [join([bit == 1 ? "X" : "I" for bit in row], ",") for row in eachrow(coordinates)]
    append!(rows, [join([bit == 1 ? "Z" : "I" for bit in row], ",") for row in eachrow(zchecks)])
    logicals = Dict{Any,Any}(
        ("T", 1) => [Op("TD", q) for q in 1:15],
        ("TD", 1) => [Op("T", q) for q in 1:15],
    )
    return StabilizerCode(rows, logicals; d=3)
end
