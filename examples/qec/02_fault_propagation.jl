# Audit every gate-local Pauli fault in a short transversal Clifford circuit.
# This uses binary Pauli propagation, so no exponentially large state is built.
# A fault occurs after one physical gate: 3 possible nonidentity Paulis for a
# one-qubit gate, or 15 for a two-qubit gate. Preparation, idle wires, readout,
# and the final ideal correction are outside this fault model.
# Run from the repository root:
#   julia --project=. examples/qec/02_fault_propagation.jl

using BlueTangle

code = steane_code(2)
decoder = steane_decoder(code)
logical_circuit = [Op("H", 1), Op("S", 2), Op("CX", 1, 2)]
physical_circuit = code.ops(logical_circuit)
paulis = ("I", "X", "Y", "Z")

# 1. Enumerate faults and propagate them through the rest of the circuit.
# Keep one result per location/fault, so failed cases remain easy to inspect.
fault_results = NamedTuple[]
for (location, gate_op) in enumerate(physical_circuit)
    targets = gate_op.q == 1 ? [gate_op.qubit] :
                              [gate_op.qubit, gate_op.target_qubit]
    for names in Iterators.product(ntuple(_ -> paulis, length(targets))...)
        all(==("I"), names) && continue
        fault = [Op(name, q; noisy=false) for (name, q) in zip(names, targets) if name != "I"]
        remaining = physical_circuit[location+1:end]
        propagated = propagate_pauli(fault, remaining; nqubits=code.n)
        block_weights = [count(op -> 7(b-1) < op.qubit <= 7b, propagated) for b in 1:code.k]
        @assert all(<=(1), block_weights)

        # Matching the syndrome alone does not establish success: a residual
        # logical Pauli also has zero syndrome. Check its logical action.
        syndrome = pauli_syndrome(code, propagated)
        recovery = decoder(syndrome)
        residual = vcat(propagated, recovery)
        logical = logical_pauli(code, residual)
        success = all(iszero, logical.x) && all(iszero, logical.z)
        push!(fault_results, (; location, gate=gate_op.name, fault=Tuple(names),
                             block_weights, logical, success))
    end
end

failures = filter(result -> !result.success, fault_results)
@assert length(fault_results) == 147
@assert isempty(failures)
println("Gate-local Pauli faults checked: ", length(fault_results))
println("Logical failures after ideal blockwise recovery: ", length(failures))

# 2. Follow one fault in the bare-ancilla syndrome circuit.
# The first check is XXXXIII, measured by an ancilla controlling four CNOTs.
# An X on that ancilla after its second CNOT propagates to the remaining data
# positions 3 and 4. The decoder adds X_7; X_3 X_4 X_7 is a logical X.
single = steane_code()
unitary_extraction = single.ops_syndrome[1:end-single.m] # omit final measurements
ancilla = single.n + 1
location = findfirst(op -> op.name == "CX" && op.qubit == ancilla &&
                          op.target_qubit == 2, unitary_extraction)
@assert location !== nothing
hook = propagate_pauli(Op("X", ancilla), unitary_extraction[location+1:end];
                       nqubits=single.n+single.m)
data_hook = filter(op -> op.qubit <= single.n, hook)
hook_syndrome = pauli_syndrome(single, data_hook)
hook_recovery = steane_decoder(single)(hook_syndrome)
hook_logical = logical_pauli(single, vcat(data_hook, hook_recovery))
@assert [(op.name, op.qubit) for op in data_hook] == [("X", 3), ("X", 4)]
@assert hook_logical == (x=[1], z=[0])
println("Bare-ancilla hook on data: ", [(op.name, op.qubit) for op in data_hook])
println("Hook syndrome: ", hook_syndrome)
println("Residual logical Pauli after ideal recovery: ", hook_logical)

# The first calculation checks this finite set of gate faults under ideal
# final recovery. The bare ancilla example fails to contain a single fault.
# Preparation, syndrome extraction, and readout require their own fault checks;
# these calculations do not estimate a threshold.
