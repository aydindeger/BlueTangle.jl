# BlueTangle.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deger.xyz/BlueTangle.jl)
[![Build Status](https://github.com/aydindeger/BlueTangle.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aydindeger/BlueTangle.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Documentation

[**BlueTangle.jl Documentation**](https://deger.xyz/BlueTangle.jl)

## Introduction

Welcome to BlueTangle.jl, a noisy quantum circuit simulator in Julia. This package is made for today's quantum computing needs and challenges. BlueTangle.jl is more than just a simulator. It's a complete toolkit for working with Noisy Intermediate-Scale Quantum (NISQ) devices.

**Key Features of BlueTangle.jl:**
- **User-Friendly and Clean Code Structure:** Thanks to Julia's multiple dispatch, BlueTangle.jl provides a straightforward and clean function namespace that's easy to use and understand.
- **Mid-Circuit Measurements:** Incorporates crucial conditional mid-circuit measurements. This feature is valuable for studying measured-based phenomena like entanglement phase transitions. It also allows for probing exotic phenomena such as the Quantum Zeno effect, where frequent measurements can effectively slow down the evolution of a quantum system.
- **Quantum Noise Channels:** Offers predefined and customizable noise channels using Kraus operators, enhancing the simulation's resemblance to real quantum environments.
- **Error Mitigation Techniques:** Features error mitigation methods, including Pauli twirling and zero-noise extrapolation, which are crucial for reducing noise in quantum computations. These techniques can be easily applied using a single function.
- **Automated Trotterized Hamiltonians:** Simplifies the creation of trotterized Hamiltonians by allowing users to input a string of operators, streamlining the process of simulating quantum systems.
- **Quantum Information Toolkit:** This toolkit enables measuring quantum entanglement, calculating partial traces, and assessing the gaussianity of states. Additionally, it features a function for calculating high-order moments, providing a unique perspective on phase transitions. It also includes simple implementation of quantum phase estimation and the ability to conduct classical shadow experiments, making it a versatile resource for quantum computing research.
- **Speed and Efficiency:** Optimized for high performance, harnessing Julia's advanced JIT compilation capabilities. This ensures fast and efficient simulations, despite being packed with a multitude of features. Its design is ideal for the iterative testing and development of complex quantum circuits, providing a balance between comprehensive functionality and swift execution.

BlueTangle.jl is designed for a diverse audience, from researchers to quantum computing enthusiasts. It stands as a robust platform for experimentation in the quantum computing landscape.

**Notice: Package Under Development**

Please be aware that this package is currently under active development and is in its alpha version. This means that the features are still being finalized, and the functionality is subject to change. Users may encounter bugs or incomplete features.

## About the Developer

Aydin Deger, a Research Fellow in Quantum Computing Theory at University College London, is the sole developer of the BlueTangle.jl. With a deep interest in theoretical physics and quantum computing, Aydin aims to bridge complex quantum concepts with practical computing applications through this package. BlueTangle.jl embodies his vision of making quantum computing more accessible and comprehensible to a diverse audience, from researchers to enthusiasts.

For further information or to explore more about Aydin's work, please visit [deger.xyz](http://deger.xyz).

## Installation
Installing BlueTangle is straightforward. You can install it through the Julia REPL using the following command:

```julia
using Pkg
Pkg.add("https://github.com/aydindeger/BlueTangle.jl")
```

Alternatively, type `] add https://github.com/aydindeger/BlueTangle.jl` directly in the Julia REPL.

To load the package in your Julia environment, use:

```julia
using BlueTangle
```

For direct installation from the source, such as in a Jupyter notebook located at the root folder of the package, use:

```julia
cd(@__DIR__)
push!(LOAD_PATH, "src")
using BlueTangle
```

## Basic usage

For more examples and advanced usage, please refer to the [**BlueTangle.jl Documentation**](https://deger.xyz/BlueTangle.jl)

### Creating and Analyzing a GHZ Circuit

This example demonstrates the creation of a Greenberger–Horne–Zeilinger (GHZ) state circuit, including the steps for measurement and calculating correlations.

```julia
# Quantum operations to prepare the GHZ state
hadamard = Op("H", 1)
cnot1 = Op("CNOT", 1, 2)
cnot2 = Op("CNOT", 2, 3)
ops = [hadamard, cnot1, cnot2]

# Compile operations into a quantum circuit
circuit = compile(ops)

# Perform measurements on the circuit
shots = 1000
measurement = sample(circuit, shots)

# Output measurement details
println("Expectation values:", measurement.expect)
println("Total magnetization moments:", measurement.mag_moments)

# Calculate and print correlations (e.g., ⟨Z₁Z₂⟩)
correlations = get_corr_from_measurement(measurement, [1, 2])
println("Correlations:", correlations)
```

## Citation 

We appreciate your support in acknowledging the use of BlueTangle.jl. For citation in research work, please use the following reference:

```
Aydin Deger, "BlueTangle.jl: A Noisy Quantum Circuit Simulator", https://github.com/aydindeger/BlueTangle/ (2024).
```

