# Functions

Detailed explanations or additional notes about all functions.

## Measurement functions

```@docs
sample_outcomes
get_probabilities_from_sample
get_expect
get_expect_from_measurement
get_corr
get_corr_from_measurement
```

## Circuit functions

```@docs
compile
quantum_circuit
sample
circuit_to_state
circuit_to_rho
```

## Plot

```@docs
plot_measurement
plot_circuit
```

## Quantum Toolkit

```@docs
entanglement_entropy
classical_shadow
mag_moments_from_rho
```

## Other

```@docs
get_N
show_basis
hilbert
hilbert_op
expand_multi_op
string_to_matrix
state_vector_create
init_state_create
apply_op!
apply_op_rho!
linear_fit
```

## Object types

```@docs
QuantumOps
Op
ifOp
Measurement
QuantumChannel
Circuit
Options
```

## Gates 

```@docs
gate
gates1
gates2
U1
U2
U3
random_ops
random_clifford
```

## Noise

```@docs
is_kraus_valid
Noise1
Noise2
apply_noise
apply_twirl
cnot_amplifier!
op_amplifier!
error_mitigate_data
```

## Hamiltonian

```@docs
trotter_ising
```

## Index

```@index
```