"""
`int2bin(a::Int, N::Int) -> Array`
Converts an integer to a binary representation as an array of bits.

These functions are utilities for bit manipulation and conversion in quantum computing contexts.
"""
# int2bin_slow(number::Int,N::Int)=reverse(digits(number,base=2,pad=N))

function int2bin(number::Int, N::Int)
    binary_vector = Vector{Int}(undef, N)
    for i in 1:N
        binary_vector[N - i + 1] = (number >> (i - 1)) & 1
    end
    return binary_vector
end


# bin2int_slow(v::Vector{Int}) = sum((2 .^ (length(v) - 1:-1:0)) .* v)
"""
`bin2int(bit_rep::Array) -> Int`
Converts a binary representation (array of bits) to an integer.
"""
function bin2int(v::Vector{Int})
    result = 0
    for bit in v
        result = (result << 1) + bit
    end
    return result
end

int_basis_create(N::Int)=collect(0:2^N-1)

fock_basis(number::Int,N::Int)=int2bin(number,N)

"""
    fock_basis_create(N::Int)
"""
fock_basis_create(N::Int)=int2bin.(0:2^N-1,N)

mag_basis(number::Int,N::Int)=1 .- 2 .* int2bin(number, N)
mag_basis_create(N::Int)=mag_basis.(0:2^N-1,N)

str2vec(str_rep::String)=[parse(Int, b, base=2) for b in str_rep]
vec2str(number_bin::Vector{Int64}) = prod(string.(number_bin))
bin2str(number_bin_vec::Vector{Vector{Int64}}) = bin2str.(number_bin_vec)

str2int(binary_string::String) = parse(Int, binary_string, base = 2)
s(number::Int,N::Int) = lpad(string(number, base=2), N, '0')

"""
`_b2v(x::Int) -> sa.SparseVector`

Converts a binary digit to a sparse vector representing quantum state.

- `x`: A binary digit (0 or 1).

Returns a sparse vector corresponding to the quantum state of the binary digit.
"""
_bin2state(x::Int)=x==0 ? sa.sparse([1.0;0.0im]) : sa.sparse([0.0im;1.0])
