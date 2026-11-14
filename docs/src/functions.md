# Functions

Reference documentation for public functions and types.

## Measurement functions

```@docs
sample
measure
born_measure_Z
get_probs_from_sample
expect
correlation
```

## Operation execution and routing

```@docs
apply
route
run(::AbstractVector{<:Number}, ::Vector{<:QuantumOps})
run(::Vector{<:QuantumOps}, ::Int)
to_state
to_rho
to_MPS
```

## Quantum Toolkit

```@docs
partial_trace
entanglement_entropy
shadow
mag_moments
```

## Plotting and optional backends

```@docs
plotq
plotqt
savefigure
phasepoly
```

## Other

```@docs
get_N
fock_basis_create
hilbert
expand_multi_op
string_to_matrix
zero_state
one_state
product_state
neel_state01
neel_state10
linear_fit
quadratic_fit
```

## Object types

```@docs
QuantumOps
Op
ifOp
Measurement
QuantumChannel
NoiseModel
Options
Layout
```

## Gates 

```@docs
gate
random_ops
random_clifford
get_stats
get_layers
```

## Noise

```@docs
is_valid_quantum_channel
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
hamiltonian
hamiltonian_exp
evolve
```

## Quantum Machine Learning

The [QML API reference](qml.md#API-reference) documents ansatz generators,
`variational_apply`, and `train` alongside executable examples.

## Index

```@index
```
