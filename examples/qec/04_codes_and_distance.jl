# Compare three codes, verify their distances, and correct a physical error.
# The code's stored d is a declaration. code_distance searches actual Pauli
# supports and supplies a witness or a proven lower bound. The correction
# below is ideal qec_round; the noisy Shor circuit currently targets Steane.
# Run from the repository root:
#   julia --project=. examples/qec/04_codes_and_distance.jl

using BlueTangle, LinearAlgebra

codes = [("Shor [[9,1,3]]", shor_code()),
         ("Steane [[7,1,3]]", steane_code()),
         ("Reed–Muller [[15,1,3]]", reed_muller_code())]

# 1. A nontrivial logical Pauli commutes with every stabilizer, yet acts on
# the encoded qubit. Searching all smaller weights certifies its minimality.
for (name, code) in codes
    distance = code_distance(code; max_weight=7)
    @assert distance.exact && distance.distance == 3
    @assert all(iszero, pauli_syndrome(code, distance.witness))
    logical_action = logical_pauli(code, distance.witness)
    @assert any(!iszero, logical_action.x) || any(!iszero, logical_action.z)
    println(name, ": d=", distance.distance, ", dX=", distance.dx,
            ", dZ=", distance.dz, ", candidates=", distance.candidates)
    println("  Minimum-weight logical witness: ",
            [(op.name, op.qubit) for op in distance.witness])
end

# 2. Stopping before finding a witness proves a bound, not the true distance.
limited = code_distance(reed_muller_code(); max_weight=2)
@assert limited.distance === nothing
@assert limited.lower_bound == 3
println("Search through weight two: d >= ", limited.lower_bound,
        ", exact distance known? ", limited.exact)

# 3. Preserve an arbitrary complex logical qubit through one physical error.
# This calculation uses ideal stabilizer projections, so it establishes code
# correction behavior rather than fault tolerance of an extraction circuit.
logical_input = normalize(ComplexF64[1, 0.3+0.6im])
for (name, code) in codes
    encoded = code.encode(logical_input)
    damaged = apply(encoded, Op("Y", 2))
    recovered = qec_round(code, damaged)
    decoded = code.decode(recovered.state)
    fidelity = abs2(dot(logical_input, decoded))
    @assert isapprox(fidelity, 1; atol=1e-10)
    println(name, " ideal corrected fidelity: ", fidelity)
end

# 4. The Reed–Muller code has an asymmetric CSS distance and a transversal T.
# Its physical TD layer implements logical T in this encoder convention.
rm = codes[3][2]
encoded = rm.encode(logical_input)
decoded_T = rm.decode(rm.apply(encoded, "T", 1))
target_T = apply(logical_input, Op("T", 1))
T_fidelity = abs2(dot(normalize(target_T), normalize(decoded_T)))
@assert isapprox(T_fidelity, 1; atol=1e-10)
println("Reed–Muller logical T fidelity: ", T_fidelity)

# Shor's nine-qubit code is a code. Shor's verified-cat extraction is a method.
# shor_qec_round(steane_code(), ...) uses that method on Steane data; it does
# not turn the data into the nine-qubit Shor code. For a noisy encoded algorithm
# and logical readout with a noisy magic-state factory, run 03_noisy_phase_estimation.jl.
