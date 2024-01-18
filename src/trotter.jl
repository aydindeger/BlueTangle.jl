"""
`trotter_ising(t::Float64, N::Int, J::Float64, h::Float64; dt=0.1) -> Vector{QuantumOps}`

Simulates the Ising model Hamiltonian using Trotterization.

- `t`: The total time for the simulation.
- `N`: Number of qubits in the Ising chain.
- `J`: Coupling constant for the interaction between adjacent qubits.
- `h`: Magnetic field strength applied in the x-direction.
- `dt` (optional): Time step for the Trotterization. Default is 0.1.

Performs Trotterization of the Ising model Hamiltonian over `N` qubits for a total time `t`, given the coupling constant `J` and magnetic field strength `h`. The function constructs a series of quantum gates that approximate the evolution under the Ising Hamiltonian.

Returns a vector of `QuantumOps`, representing the sequence of operations for the Trotterized Ising model simulation.
"""
function trotter_ising(t::Float64,N::Int,J::Float64,h::Float64;dt=0.1)
    trotter_step=round(Int,t/dt)
    x_angle=-2*dt*h
    z_angle=-2*dt*J
    final_list=vcat([Op("RX($(x_angle))",i) for i=1:N],[[Op("CNOT",i,i+1),Op("RZ($(z_angle))",i+1),Op("CNOT",i,i+1)] for i=1:N-1]...)
    
    return vcat(fill(final_list,trotter_step)...)#list of operators
end