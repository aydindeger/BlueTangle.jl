# A fault-tolerant logical iSWAP between two [[6,2,2]] blocks, compared with bare qubits and with one
# [[4,2,2]] block. Run from the repository root:
#   julia --project=. examples/qec/07_fault_tolerant_iswap.jl
#
# 1. The circuit: prepare logical inputs in two blocks, apply k logical iSWAPs, read every wire out in one basis per
#    block. BlueTangle's statevector simulator confirms the predicted logical outcomes; every shot passes the checks.
# 2. Every single fault, one at a time: none is harmful (undetected and changing a logical bit).
# 3. The logical error after 8 logical iSWAPs against the physical error p of their gates, for three schemes
#    (preparation and readout ideal, every gate of the iSWAPs noisy):
#      bare qubits:  physical iSWAPs                            -> error ∝ p
#      [[4,2,2]]:    physical iSWAPs inside one block            -> undetected error ∝ p (faults XX, YY, ZZ)
#      [[6,2,2]]:    transversal iSWAPs between two blocks        -> undetected error ∝ p²
# 4. The whole circuit with noise everywhere, as an experiment would see it.

using BlueTangle, LinearAlgebra, Printf

# 1. The circuit and its ideal outcomes.
k = 3
flip = Dict(:Z => :Z, :X => :Y, :Y => :X)
inputs_a, inputs_b = (1, 0), (0, 1)         # logical values of A1, A2 (X basis) and B1, B2 (Z basis)
out = isodd(k) ? (flip[:Z], flip[:X]) : (:X, :Z)                     # (Z, Y) for odd k
circuit = vcat(c6_prepare(:X; bits=inputs_a, block=1), c6_prepare(:Z; bits=inputs_b, block=2),
               reduce(vcat, [c6_iswap() for _ in 1:k]), c6_readout(out[1]; block=1), c6_readout(out[2]; block=2))
psi = apply(zero_state(12), circuit)
results = [c6_decode(int2bin(s, 12), out) for s in sample(psi, 500)]
@assert all(r.accepted for r in results)
@assert length(unique(r.logical for r in results)) == 1
println("1. ", k, " logical iSWAPs on |-,+> (A, X basis) and |0,1> (B, Z basis): every shot accepted; logical ",
        "outcome A in ", out[1], ", B in ", out[2], ": ", first(results).logical)
println("   Gates: ", count(o -> o.name == "CX", circuit), " CNOTs in the two encoders, ",
        count(o -> o.name == "ISWAP", circuit), " physical iSWAPs (two native sqrt-iSWAPs each)")

# 2. Every single fault.
audit = c6_single_faults(circuit)
@assert audit.harmful == 0
println("2. ", audit.faults, " single faults: ", audit.detected, " detected, ", audit.harmless, " harmless, ",
        audit.harmful, " harmful")

# 3. Eight logical iSWAPs, noisy gates, ideal preparation and readout.
#    The final Pauli frame decides: detected (some check anticommutes), or a logical error (undetected and not a
#    stabilizer), counted per logical iSWAP pair, so that the three schemes are compared on the same footing.
const C6_CHECKS = ([1, 0, 1, 1, 1, 0], [0, 1, 1, 0, 1, 1])
const C6_LOGICAL = ([1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1])
parity(bits, support) = isodd(sum(bits .* support))

function c6_pair_errors(x, z)
    # Per shot: detected? and, if not, which logical pairs (A_i, B_i) carry a logical error.
    detected = falses(size(x, 1))
    wrong = falses(size(x, 1), 2)
    for s in axes(x, 1), b in 1:2
        xs, zs = x[s, 6(b-1)+1:6b], z[s, 6(b-1)+1:6b]
        detected[s] |= any(parity(xs, c) || parity(zs, c) for c in C6_CHECKS)
        for i in 1:2
            wrong[s, i] |= parity(xs, C6_LOGICAL[i]) || parity(zs, C6_LOGICAL[i])
        end
    end
    kept = .!detected
    return sum(wrong[kept, :]) / (2 * max(sum(kept), 1))
end

function c4_errors(x, z)
    # One [[4,2,2]] block: checks XXXX and ZZZZ detect an odd X or Z weight; an undetected frame is harmless only
    # if it is a stabilizer (X part and Z part each 0000 or 1111).
    detected = isodd.(sum(x; dims=2)[:]) .| isodd.(sum(z; dims=2)[:])
    trivial(v) = all(v) || !any(v)
    wrong = [!(trivial(x[s, :]) && trivial(z[s, :])) for s in axes(x, 1)]
    kept = .!detected
    return sum(wrong[kept]) / max(sum(kept), 1)
end

bare_errors(x, z) = sum(any(x .| z; dims=2)) / size(x, 1)

K = 8
gadget_c6 = reduce(vcat, [[Op("ISWAP", j, j + 6) for j in 1:6] for _ in 1:K])   # transversal; Z layers are free
gadget_c4 = [Op("ISWAP", 2, 3) for _ in 1:K]                                     # inside one [[4,2,2]] block
gadget_bare = [Op("ISWAP", 1, 2) for _ in 1:K]
ps = [1e-3, 2e-3, 4e-3, 8e-3, 1.6e-2]
table = Dict(:bare => Float64[], :c4 => Float64[], :c6 => Float64[])
shots = 1_000_000
for p in ps
    f = sample_pauli_frames(gadget_bare, p; nqubits=2, shots=shots)
    push!(table[:bare], bare_errors(f.x, f.z))
    f = sample_pauli_frames(gadget_c4, p; nqubits=4, shots=shots)
    push!(table[:c4], c4_errors(f.x, f.z))
    f = sample_pauli_frames(gadget_c6, p; nqubits=12, shots=shots)
    push!(table[:c6], c6_pair_errors(f.x, f.z))
end
function slope(v)
    good = v .> 0
    A = hcat(ones(count(good)), log.(ps[good]))
    return (A \ log.(v[good]))[2]
end
println("3. $K logical iSWAPs, logical error per pair, every gate noisy with probability p ($shots shots per point):")
println("        p      bare     [[4,2,2]]   [[6,2,2]] x 2")
for (i, p) in enumerate(ps)
    @printf("   %7.0e   %.2e   %.2e    %.2e\n", p, table[:bare][i], table[:c4][i], table[:c6][i])
end
slopes = Dict(key => slope(v) for (key, v) in table)
@printf("   log-log slopes: bare %.2f, [[4,2,2]] %.2f, [[6,2,2]] %.2f\n", slopes[:bare], slopes[:c4], slopes[:c6])
@assert isapprox(slopes[:bare], 1; atol=0.1) && isapprox(slopes[:c4], 1; atol=0.15) && isapprox(slopes[:c6], 2; atol=0.3)
println("   Only the transversal [[6,2,2]] iSWAP is fault tolerant: its undetected error is second order in p.")

# 4. The whole circuit, every gate noisy (error p per physical iSWAP, p/2 per one-qubit gate, readout error p), as
#    an experiment on the 2 x 6 ladder would see it: kept fraction and logical error of kept shots against k, next
#    to bare qubits running the same iSWAP layers.
p = 3e-3
rotation(basis, qubits) = basis == :Z ? Op[] : basis == :X ? Op[Op("H", q) for q in qubits] :
                          vcat(Op[Op("SD", q) for q in qubits], Op[Op("H", q) for q in qubits])
println("4. Whole circuit, p = $p per iSWAP, p/2 per single-qubit gate, readout p:")
println("    k   kept    encoded error   bare error")
for k in (0, 1, 2, 4, 8, 16)
    bases = isodd(k) ? (:Z, :Y) : (:X, :Z)               # inputs: block A in X, block B in Z
    layers = reduce(vcat, [c6_iswap() for _ in 1:k]; init=Op[])
    enc = vcat(c6_prepare(:X; block=1), c6_prepare(:Z; block=2), layers,
               c6_readout(bases[1]; block=1), c6_readout(bases[2]; block=2))
    f = sample_pauli_frames(enc, p / 2; nqubits=12, shots=200_000, p2=p, readout=p)
    kept, wrong = 0, 0
    for s in axes(f.x, 1)
        r = c6_decode(Int.(f.x[s, :]), (:Z, :Z))          # parities of the flips (the basis only adds an offset)
        r.accepted || continue
        kept += 1
        wrong += sum(r.logical)
    end
    bare = vcat(Op[Op("H", q) for q in 1:6], reduce(vcat, [[Op("ISWAP", j, j + 6) for j in 1:6] for _ in 1:k]; init=Op[]),
                rotation(bases[1], 1:6), rotation(bases[2], 7:12))
    g = sample_pauli_frames(bare, p / 2; nqubits=12, shots=200_000, p2=p, readout=p)
    @printf("   %2d   %.2f      %.2e        %.2e\n", k, kept / 200_000, wrong / (4kept), sum(g.x) / (12 * 200_000))
end
println("Encoded errors of kept shots stay well below the bare errors; the kept fraction falls with k.")
