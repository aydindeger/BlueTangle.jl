# Prepare an encoded Bell pair, propagate a physical fault through a logical
# CNOT, and recover the two blocks. Encoding and correction are ideal here;
# only the transversal logical gates have the demonstrated fault-containment
# property. This script performs local simulation and submits no hardware jobs.
# Run from the repository root:
#   julia --project=. examples/qec/01_logical_bell.jl

using BlueTangle, LinearAlgebra

# Each logical qubit occupies seven consecutive physical qubits.
code = steane_code(2)
decoder = steane_decoder(code)
logical_circuit = [Op("H", 1), Op("CX", 1, 2)]
fault = Op("X", 2)                 # physical qubit 2, in the control block

# 1. Establish the ideal logical answer and its encoded implementation.
logical_input = zero_state(2)
expected = apply(logical_input, logical_circuit)
encoded_input = code.encode(logical_input)
physical_circuit = code.ops(logical_circuit)
ideal_encoded = apply(encoded_input, physical_circuit)
@assert isapprox(abs2(dot(expected, code.decode(ideal_encoded))), 1; atol=1e-12)

# 2. Insert X before the transversal CNOT. CNOT maps X on its control to
# X on both operands: the output error is X_2 X_9. Each block still contains
# only one erroneous qubit, although the complete register has two errors.
prepared = apply(encoded_input, code.ops("H", 1))
damaged = apply(prepared, fault)
damaged = apply(damaged, code.ops("CX", 1, 2))
output_fault = propagate_pauli(fault, code.ops("CX", 1, 2); nqubits=code.n)
@assert [(op.name, op.qubit) for op in output_fault] == [("X", 2), ("X", 9)]
println("Output physical fault: ", [(op.name, op.qubit) for op in output_fault])

# 3. Recover each Steane block. qec_round measures ideal stabilizer projectors
# and returns only the data qubits, so correction rounds can be composed.
# A decoder that searches for one error in the entire register is insufficient
# here; steane_decoder treats each seven-qubit block independently.
first_round = qec_round(code, damaged; decoder=decoder)
decoded = code.decode(first_round.state)
bell_fidelity = abs2(dot(expected, decoded))
@assert get_N(first_round.state) == 14
@assert isapprox(bell_fidelity, 1; atol=1e-12)
println("First syndrome: ", first_round.syndrome)
println("Recovered logical Bell fidelity: ", bell_fidelity)

# 4. Continue the computation, inject another error, and correct again.
# These extra gates also check a complex logical phase: on this Steane code,
# logical S is implemented by physical SD on all seven qubits of its block.
next_gates = [Op("S", 1), Op("H", 2)]
expected_next = apply(expected, next_gates)
continued = apply(first_round.state, code.ops(next_gates))
continued = apply(continued, Op("Z", 10))
second_round = qec_round(code, continued; decoder=decoder)
continued_fidelity = abs2(dot(expected_next, code.decode(second_round.state)))
@assert get_N(second_round.state) == 14
@assert isapprox(continued_fidelity, 1; atol=1e-12)
println("Second syndrome: ", second_round.syndrome)
println("Recovered continued-circuit fidelity: ", continued_fidelity)
println("Data qubits after both rounds: ", get_N(second_round.state))

# Change logical_circuit/next_gates to build a Clifford algorithm. Edit fault
# locations to investigate propagation. For a circuit-level fault-tolerance
# claim, replace ideal preparation and qec_round with validated noisy gadgets;
# the second example demonstrates why the existing bare-ancilla circuit fails.
