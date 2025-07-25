module BlueTangle

include("all.jl")

export get_N, ro3, ro10, fields, sample, sample_bit, get_probs_from_sample, sample_state, sample_exact, expand_multi_op, string_to_matrix
export expect, correlation, apply, measure, compile, quantum_circuit, prob, to_state, to_rho, fock_basis
export QuantumOps, Op, ifOp, OpF, OpQC, Measurement, QuantumChannel, Layout, Circuit, Options, NoiseModel #Structs
export init, N_MPS, hilbert, hilbert_op, product_state, zero_state, one_state, plus_state, minus_state, neel_state01, neel_state10, random_product_state, random_state, neel_state_s, Z3, Z3_s
export pauli_decomposition,pauli_reconstruction, pauli_decomposition_names, pauli_decomposition_tensor
export gate,gates,random_ops,random_ops_len,random_clifford,Noise1,Noise2,apply_noise,U1,U2,U3,is_valid_quantum_channel,apply_twirl,custom_noise,cnot_amplifier!,op_amplifier!,linear_fit,quadratic_fit,error_mitigate_data
export plotq, savefigure #Plot
export ⊗, isunitary, same_up_to_global_phase, partial_trace, stinespring_dilation, reduced_row_echelon, rank_of_rref #Linear algebra, linalg
export fock_basis_create, int2bin, bin2int, sparsevector, hamming_distance, born_measure_Z
export zyz_decomposition, kronecker_decomposition #Decomposition
export maxdim, _born_measure
export StabilizerCode, encoding_circuit_from_generator, physical_ops, qec_state_prep, QECState, get_XZ_logicals!, get_XZ_logicals, stabilizers_to_generator, reduce_and_relabel_qubits, get_standard_form #QEC get_codeword, 
export relabel_swap
export entanglement_entropy, shadow, mag_moments
export to_MPS, to_state, inner, amplitude
export it, its, sa, plt, la, sb #Packages
export hamiltonian, evolve, hamiltonian_exp, AnsatzOptions, VQE, variational_apply, get_stats, get_layers, measure_ZNE, fidelity, get_optimized_layers, optimize_simple

end