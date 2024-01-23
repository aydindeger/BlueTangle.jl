"""
`int2bit(a::Int, N::Int) -> Array`
Converts an integer to a binary representation as an array of bits.

These functions are utilities for bit manipulation and conversion in quantum computing contexts.
"""
int2bit(a::Int,N::Int)=reverse(digits(a,base=2,pad=N))

"""
`show_basis(N::Int)=int2bit.(0:2^N-1,N)`
"""
show_basis(N::Int)=int2bit.(0:2^N-1,N)

"""
`bit2label(bit_rep::Array) -> Int`
Converts a binary representation (array of bits) to an integer.
"""
function bit2label(bit_rep)
    digit=0
    for (i,bit) in enumerate(reverse(bit_rep))
        digit += bit*2^(i-1)
    end
    return Int(digit) #state_number can be 0
  end

"""
`str2label(str_rep::String) -> Int`
Converts a binary representation (string) to an integer.
"""
str2label(str_rep)=[parse(Int, b, base=2) for b in str_rep]


#todo make this more efficient
"""
`bit_to(bit_basis::Vector, N::Int, sym::Symbol) -> Any`

Converts bit basis representations to different formats based on the symbol provided.

- `bit_basis`: Vector representing the bit basis.
- `N`: Number of bits.
- `sym`: Symbol indicating the desired format (:fock, :mag, :bstr, etc.).

Returns the bit basis in the specified format.
"""
function bit_to(bit_basis::Vector,N::Int,sym::Symbol)
    # bit_basis=0:2^N-1
    # N=Int(log(2,length(bit_basis)))
    d=Dict()
    fock=int2bit.(bit_basis,N)
    d[:int]=bit_basis
    d[:fock]=fock
    d[:mag]=replace.(fock,1=>-1,0=>1)
    d[:bstr]=join.(map.(string, fock))

    return d[sym]
end

bit_to(bit_basis::UnitRange,N::Int,sym::Symbol)=bit_to(collect(bit_basis),N,sym)
bit_to(N::Int,sym::Symbol)=bit_to(collect(0:2^N-1),N,sym)

"""
`_b2v(x::Int) -> SparseVector`

Converts a binary digit to a sparse vector representing quantum state.

- `x`: A binary digit (0 or 1).

Returns a sparse vector corresponding to the quantum state of the binary digit.
"""
function _b2v(x)
    if x==0 || x==1
        return x==0 ? sparse([1.0;0.0im]) : sparse([0.0im;1.0])
    elseif x==[1;0] || x==[0;1]
        return x==[1;0] ? 0 : 1
    end
end

#todo make this more efficient
"""
`state_vector_create(list_of_qubits::Vector) -> SparseVector`

Creates a quantum state vector from a list of qubits.

- `list_of_qubits`: A vector representing the state of each qubit.

Returns a sparse vector representing the quantum state of the system.
"""
state_vector_create(list_of_qubits::Vector)=sparse(foldl(kron,_b2v.(list_of_qubits)))

"""
`init_state_create(N::Int) -> SparseVector`

Returns a sparse vector representing the |000...> quantum state.
"""
init_state_create(N::Int)=spzeros(ComplexF64,2^N);state[1]=1;