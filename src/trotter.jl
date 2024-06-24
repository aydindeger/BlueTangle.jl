"""
    trotter_ising(t::Float64, N::Int, J::Float64, h::Float64; dt=0.1) -> Vector{QuantumOps}

Simulates the Ising model Hamiltonian using Trotterization.

- `total_time`: The total time for the simulation.
- `N`: Number of qubits in the Ising chain.
- `J`: Coupling constant for the interaction between adjacent qubits.
- `h`: Magnetic field strength applied in the x-direction.
- `dt` (optional): Time step for the Trotterization. Default is 0.1.

Performs Trotterization of the Ising model Hamiltonian over `N` qubits for a total time `total_time`, given the coupling constant `J` and magnetic field strength `h`. The function constructs a series of quantum gates that approximate the evolution under the Ising Hamiltonian.

Returns a vector of `QuantumOps`, representing the sequence of operations for the Trotterized Ising model simulation.
"""
function trotter_ising(N::Int,total_time::Float64,J::Float64,h::Float64;dt=0.1)
    trotter_step=round(Int,total_time/dt)
    x_angle=round(-2*dt*h,sigdigits=6)
    z_angle=round(-2*dt*J,sigdigits=6)
    final_list=vcat([Op("RX($(x_angle))",i) for i=1:N],[[Op("CNOT",i,i+1),Op("RZ($(z_angle))",i+1),Op("CNOT",i,i+1)] for i=1:N-1]...)
    
    return vcat(fill(final_list,trotter_step)...)#list of operators
end

"""
    hamiltonian_exp(N::Int, total_time::Float64, string_of_ops::Vector; dt=0.1) -> Vector{QuantumOps}

Expands a given Hamiltonian expressed as a string of operations into a sequence of quantum gates using Trotterization.

- `N`: Number of qubits involved in the simulation.
- `total_time`: The total simulation time over which the Hamiltonian is to be applied.
- `string_of_ops`: A vector where odd indices contain the coupling strengths and even indices contain comma-separated strings representing the operators (X, Y, Z) applied to consecutive qubits.
- `dt` (optional): The time step for Trotterization, with a default value of 0.1.

The function `hamiltonian_exp` parses the `string_of_ops` to construct a sequence of operations based on the specified operators and their coupling strengths. For each term in the Hamiltonian:

This method returns a vector of `QuantumOps`, each representing a quantum operation to be applied sequentially to simulate the Hamiltonian over the specified time `total_time`.

### Example
# Define a Hamiltonian for a 3-qubit system with mixed interactions over 2 seconds
N = 3
total_time = 2.0
string_of_ops = [1.0, "X,Y", 0.5, "Y,Z"]
ops = hamiltonian_exp(N, total_time, string_of_ops)

# ops will contain a sequence of quantum gates to apply.
"""
function hamiltonian_exp(N::Int,total_time::Float64,string_of_ops::Vector;dt=0.01,full=true)

    couplings=[]
    ops=[]
    for (i,op)=enumerate(string_of_ops)
        if isodd(i)
            push!(couplings,op)
        else
            push!(ops,op)
        end
    end
    
    term_ops=[split(op,",") for op=ops]
    len_op=[length(split(op,",")) for op=ops]
    
    trotter_step=round(Int,total_time/dt)
    
    all_ops=Vector{QuantumOps}()
    
    for term=1:length(ops)
    
        terms=term_ops[term]
        j=len_op[term]-1
        angle=round(-2*dt*couplings[term],sigdigits=6)
    
        for i=1:N-j
    
            for (en,k)=enumerate(i:i+j)
                if terms[en]=="X"
    
                    if j==0
                        push!(all_ops,Op("RX($(angle))",i+j))
                    else
                        push!(all_ops,Op("H",k))
                    end
    
                elseif terms[en]=="Y"
    
                    if j==0
                        push!(all_ops,Op("RY($(angle))",i+j))
                    else
                        # push!(all_ops,Op("HSP",gate.HSP,k)) 
                        push!(all_ops,Op("HY",k)) 
                    end
    
                elseif terms[en]=="Z"
    
                    if j==0
                        push!(all_ops,Op("RZ($(angle))",i+j))
                    end

                else

                    throw("Only X,Y,Z operations are allowed")
                    
                end
            end
    
            for k=i:i+j-1
                push!(all_ops,Op("CNOT",k,k+1))
            end
    
            if j>0
                push!(all_ops,Op("RZ($(angle))",i+j))
            end
    
            for k=i+j-1:-1:i
                push!(all_ops,Op("CNOT",k,k+1))
            end
    
            for (en,k)=enumerate(i:i+j)
                if terms[en]=="X" && j>0
                    push!(all_ops,Op("H",k))
                elseif terms[en]=="Y" && j>0
                    # push!(all_ops,Op("HSP'",Matrix(gate.HSP'),k))
                    push!(all_ops,Op("HY",k))
                end
            end
    
        end
    
    end
    
    return full ? vcat(fill(all_ops,trotter_step)...) : all_ops
    
    end

##

