"""
N_MPS(N::Int)
"""
N_MPS(N::Int)=it.siteinds("Qubit", N)

##
##========== state preparation ==========

"""
create all zero state
"""
zero_state(M::Vector)=it.productMPS(M,"0")

"""
create all one state
"""
one_state(M::Vector)=it.productMPS(M,"1")

"""
create given product state
"""
product_state(M::Vector,list_of_qubits::Vector)=it.productMPS(M,map(string,Int.(list_of_qubits)))

"""
create neel state 010101
"""
neel_state01(M::Vector)=product_state(M,[isodd(i) ? 0 : 1 for i=1:length(M)])

"""
create neel state 101010
"""
neel_state10(M::Vector)=product_state(M,[isodd(i) ? 1 : 0 for i=1:length(M)])

"""
create random state
"""
random_state(sites::Vector,chi::Int=1)=it.randomMPS(sites,chi)


"""
creates MPS plus state
"""
plus_state(n::Vector)=convert_X_basis(zero_state(n))

"""
creates MPS minus state
"""
minus_state(n::Vector)=convert_X_basis(one_state(n))

function convert_X_basis(state::Union{AbstractVectorS,it.MPS})
    N=get_N(state)
    for i=1:N
        state=Op("H",i)*state
    end
    return state
end

##========== state preparation ==========

maxlinkdim(MPS::it.MPS)=it.maxlinkdim(MPS)
maxlinkdim(MPO::it.MPO)=it.maxlinkdim(MPO)


expect(psi::it.MPS, op_str::String)=it.expect(psi, op_str)

expect(psi::it.MPS, mpo::it.MPO)=it.inner(psi', mpo, psi)


# ψ0 = MPS(s, n -> n == 1 ? "↓" : "↑")

"""
to_state(MPS::it.MPS,M::Vector)

    Convert MPS to state vector
"""
to_state(MPS::it.MPS,M::Vector)=sa.sparse(reshape(its.array(its.contract(MPS),reverse(M)),2^length(M)))
to_state(MPS::it.MPS)=sa.sparse(reshape(its.array(its.contract(MPS),reverse(its.siteinds(MPS))),2^length(MPS)))

function reflectMPS(psi::it.MPS)
    N = length(psi)
    psi_reflected = it.MPS(N)
    for i in 1:N
        A = psi[N - i + 1]
        psi_reflected[i] = A#it.dag(A)
    end
    return psi_reflected
end

"""
    to_MPS(vec::AbstractVector, sites::AbstractVector; kwargs...)

    Convert state vector to MPS
"""
function to_MPS(vec::AbstractVectorS, sites::AbstractVector; kwargs...)
    length(vec) == 2^length(sites) || throw(DimensionMismatch("state size must match the qubit sites"))
    # Julia's first tensor index runs fastest; qubit 1 is our most significant bit.
    tensor = its.ITensor(reshape(Vector(vec), ntuple(_ -> 2, length(sites))), reverse(sites)...)
    return it.MPS(tensor, sites; kwargs...)
end

function tensor_to_matrix(tensor::its.ITensor)
    a=tensor.tensor.storage
    len=length(a)
    if len==2
        return a #pure
    else
        s=Int(log2(len))
        return reshape(a,(s,s))
    end
end



"""
    amplitude(psi::it.MPS,config::Vector)

this will give amplitude of given configuration

Example:
N=3
config=[0,0,0] or [1,1,1]
"""
function amplitude(psi::it.MPS,config::Vector)

    M=it.siteinds(psi)
    N=length(M)

    if length(config) != N
        throw("config must be same size as system size")
    end

    config=config .+1
  
    V = it.ITensor(1.)
    for j=1:N
      V *= (psi[j]*it.state(M[j],config[j]))
    end
  
    return it.scalar(V)

end

"""
    prob(state::Union{sa.SparseVector,it.MPS},shots::Int,select_qubit::Int)
this gives probability of measuring 0 and 1 on select_qubit after shots
"""
function prob(state::Union{AbstractVectorS,it.MPS},shots::Int,select_qubit::Int)

    shots > 0 || throw(ArgumentError("shots must be positive"))
    1 <= select_qubit <= get_N(state) || throw(ArgumentError("invalid selected qubit"))

    all_sample=hcat(sample_bit(state,shots)...)[select_qubit,:]

    P1=sum(all_sample)/shots
    P0=1-P1

    return P0,P1
end


"""
sample state/MPS
"""
function sample_bit(state::AbstractVectorS,shots::Int=1;basis="Z")
    shots >= 0 || throw(ArgumentError("shots must be nonnegative"))
    if uppercase(basis)=="X"
        state=convert_X_basis(state)
    elseif uppercase(basis)=="Z"
        nothing
    else
        throw("basis must be Z or X")
    end

    return int2bin.(sample(state,shots),get_N(state))

end


function sample_bit(psi::it.MPS,shots::Int=1;basis="Z")

    shots >= 0 || throw(ArgumentError("shots must be nonnegative"))
    MPS = copy(psi)

    if uppercase(basis)=="X"
        MPS=convert_X_basis(MPS)
    elseif uppercase(basis)=="Z"
        nothing
    else
        throw("basis must be Z or X")
    end
    
    # Prepare the sampling gauge once, and leave the caller's MPS untouched.
    it.orthogonalize!(MPS, 1)
    return [it.sample(MPS).-1 for i=1:shots]
end


"""
inner(ψ::AbstractVectorS,MPS::it.MPS)=ψ'to_state(MPS)
inner(MPS::it.MPS,ψ::AbstractVectorS)=conj(inner(ψ,MPS))
inner(ψ::AbstractVectorS,ψ2::AbstractVectorS)=ψ'ψ2
inner(MPS::it.MPS,MPS2::it.MPS)=it.inner(MPS,MPS2)
"""
inner(ψ::AbstractVectorS,MPS::it.MPS)=ψ'to_state(MPS)
inner(MPS::it.MPS,ψ::AbstractVectorS)=conj(inner(ψ,MPS))
inner(ψ::AbstractVectorS,ψ2::AbstractVectorS)=ψ'ψ2
inner(MPS::it.MPS,MPS2::it.MPS)=it.inner(MPS,MPS2)


function _mat_to_tensor(sites::Vector,mat::AbstractMatrix,qubit::Int,target_qubit::Int;control::Int=-2)

    if control==-2
        return it.op(mat,sites[target_qubit],sites[qubit])#note how target qubit comes first. this is correct!#note how target qubit comes first. this is correct!
    else
        return it.op(mat,sites[target_qubit],sites[control],sites[qubit])
    end
end

_mat_to_tensor(sites::Vector,mat::AbstractMatrix,qubit::Int)=it.op(mat,sites[qubit])

"""
    fidelity(ψ::AbstractVectorS,ψ2::AbstractVectorS)
"""
fidelity(ψ::AbstractVectorS,ψ2::AbstractVectorS)=abs2(ψ'ψ2)


fidelity(rho::sa.SparseMatrixCSC, sigma::sa.SparseMatrixCSC)=fidelity(Matrix(rho),Matrix(sigma))
function fidelity(rho::AbstractMatrix, sigma::AbstractMatrix)
    a = sqrt(rho) * sigma * sqrt(rho)
    b = sqrt(a)
    fidelity_value = la.tr(b)^2
    
    return real(fidelity_value)
end

"""
    fidelity(ψ::it.MPS,ψ::AbstractVectorS)=abs2(inner(ψ,MPS))
"""
fidelity(ψ::it.MPS,ψ2::it.MPS)=abs2(inner(ψ,ψ2))

function entanglement_entropy(psi::it.MPS, b::Int; spectrum_bool::Bool=false)
    ψ = copy(psi)
    N = length(ψ)
    @assert 1 ≤ b < N "b must satisfy 1 ≤ b < N (cut between b and b+1)."

    s = it.siteinds(ψ)
    it.orthogonalize!(ψ, b)

    if b == 1
        left_inds = (s[b],)  # no left link at site 1
    else
        left_inds = (it.linkind(ψ, b-1), s[b])
    end
    _, S, _ = it.svd(ψ[b], left_inds)

    d = it.dim(it.commonind(S))
    p = Vector{Float64}(undef, d)
    for n in 1:d
        p[n] = abs2(S[n,n])          # <-- safer than (S[n,n])^2
    end

    # normalise (crucial!)
    Z = sum(p)
    if Z == 0.0
        return spectrum_bool ? (0.0, Float64[]) : 0.0
    end
    p ./= Z

    # remove/clamp numerical junk (more robust than just > 0.0)
    p = p[p .> 0]

    SvN = -sum(p .* log.(p))

    if spectrum_bool
        spectrum = sort(-log.(p))    # sort helps many downstream routines
        return SvN, spectrum
    else
        return SvN
    end
end

function entanglement_entropy(psi::it.MPS; spectrum_bool::Bool=false)
    length(psi) == 1 && return spectrum_bool ? (0.0, [0.0]) : 0.0
    return entanglement_entropy(psi, length(psi) ÷ 2; spectrum_bool=spectrum_bool)
end

function expect(psi::it.MPS, op::Op)
    applied = it.apply(op.expand(get_M(psi)), psi; cutoff=0.0)
    return real(it.inner(psi, applied))
end
expect(psi::it.MPS, name::String, qubit::Int) = expect(psi, Op(name, qubit))

function correlation(psi::it.MPS, names::String, qubits::Vector)
    operators = strip.(split(names, ','))
    length(operators) == length(qubits) || throw(DimensionMismatch("one operator is required per qubit"))
    result = copy(psi)
    sites = get_M(psi)
    for (name, qubit) in reverse(collect(zip(operators, qubits)))
        result = it.apply(it.op(gates(String(name)), sites[qubit]), result; cutoff=0.0)
    end
    return real(it.inner(psi, result))
end
correlation(psi::it.MPS, qubits::Vector) = correlation(psi, join(fill("Z", length(qubits)), ','), qubits)
