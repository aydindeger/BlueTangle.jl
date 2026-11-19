# Teleport an encoded qubit with a physical Y error on its source block.
# The two Bell readouts locate the error and preserve the logical amplitudes.
# Run locally with julia --project=. examples/qec/07_logical_teleportation.jl

using BlueTangle, LinearAlgebra

code = steane_code()
logical = normalize(ComplexF64[1, 0.2 + 0.7im])
source = code.encode(logical)       # known input for the fidelity check
fault_site = 3                      # choose any physical wire from 1 to 7
damaged = apply(source, Op("Y", fault_site))

# The gadget prepares both halves of the logical Bell pair. A source Y has
# a Z component that flips one X readout and an X component that propagates
# through the source-to-sender CNOT and flips one sender Z readout.
result = steane_teleport(code, damaged)
for (basis, readout) in (("source X", result.x_readout), ("sender Z", result.z_readout))
    println(basis, " raw bits = ", readout.bits, ", corrected wire = ",
            readout.corrected_qubit, ", logical outcome = ", readout.outcome)
    @assert readout.corrected_qubit == fault_site
end

# Both measured blocks have been consumed. Only the seven receiver wires
# remain. Ideal decoding here lets us compare to the original complex state.
fidelity = abs2(dot(logical, code.decode(result.state)))
@assert get_N(result.state) == 7
@assert isapprox(fidelity, 1; atol=1e-9)
println("Teleported logical-state fidelity = ", fidelity)
