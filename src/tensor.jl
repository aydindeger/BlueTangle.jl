"""
N_MPS(N::Int)
"""
N_MPS(N::Int)=it.siteinds("Qubit", N)

"""
init(N::Int)
"""
init(N::Int)=N,sitesN(N)

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

##========== state preparation ==========

dim(MPS::it.MPS)=it.maxlinkdim(MPS)
dim(MPO::it.MPO)=it.maxlinkdim(MPO)


expect(psi::it.MPS,op_str::String)=it.expect(psi, op_str)


# ψ0 = MPS(s, n -> n == 1 ? "↓" : "↑")

"""
to_state(MPS::it.MPS,M::Vector)

    Convert MPS to state vector
"""
to_state(MPS::it.MPS,M::Vector)=sa.sparse(reshape(its.array(its.contract(MPS),reverse(M)),2^length(M)))
to_state(MPS::it.MPS)=sa.sparse(reshape(its.array(its.contract(MPS),reverse(its.siteinds(MPS))),2^length(MPS)))

#this is wrong!
# to_state(MPS::it.MPS)=round.(sa.sparse(reshape(it.contract(MPS).tensor,2^length(MPS))),digits=10)

function reflectMPS(psi::it.MPS)
    N = length(psi)
    psi_reflected = it.MPS(N)
    for i in 1:N
        A = psi[N - i + 1]
        psi_reflected[i] = A#it.dag(A)
    end
    return psi_reflected
end

# to_MPS(vec::Vector,M::Vector,maxdim::Int=500)=it.MPS(vec,M;maxdim=maxdim,cutoff=1e-10)|>reflectMPS
to_MPS(vec::Vector,M::Vector)=it.MPS(vec,M)|>reflectMPS

# to_MPS(vec::AbstractVectorS,M::Vector,maxdim::Int=500)=to_MPS(Vector(vec),M,maxdim)

"""
    to_MPS(vec::AbstractVectorS,M::Vector)

    Convert state vector to MPS
"""
to_MPS(vec::AbstractVectorS,M::Vector)=to_MPS(Vector(vec),M)

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
function prob(state::Union{sa.SparseVector,it.MPS},shots::Int,select_qubit::Int)

    all_sample=hcat(sample_bit(state,shots)...)[select_qubit,:]

    P1=sum(all_sample)/shots
    P0=1-P1

    return P0,P1
end


"""
sample state/MPS
"""
sample_bit(state::AbstractVectorS,shots::Int=1)=int2bin.(sample(state,shots),get_N(state))
sample_bit(MPS::it.MPS,shots::Int=1)=[it.sample!(MPS).-1 for i=1:shots]

"""
inner(ψ::AbstractVectorS,MPS::it.MPS)=ψ'to_state(MPS)
inner(MPS::it.MPS,ψ::AbstractVectorS)=inner(ψ,MPS)
inner(ψ::AbstractVectorS,ψ2::AbstractVectorS)=ψ'ψ2
inner(MPS::it.MPS,MPS2::it.MPS)=it.inner(MPS',MPS2)
"""
inner(ψ::AbstractVectorS,MPS::it.MPS)=ψ'to_state(MPS)
inner(MPS::it.MPS,ψ::AbstractVectorS)=inner(ψ,MPS)
inner(ψ::AbstractVectorS,ψ2::AbstractVectorS)=ψ'ψ2
inner(MPS::it.MPS,MPS2::it.MPS)=it.inner(MPS',MPS2)


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

# inner_slow(MPS::it.MPS,ψ::AbstractVectorS,maxdim::Int)=it.inner(MPS',to_MPS(ψ,it.siteinds(MPS),maxdim))
# inner_slow(MPS::it.MPS,ψ::AbstractVectorS)=it.inner(MPS',to_MPS(ψ,it.siteinds(MPS)))

# function entanglement_entropy(psi::it.MPS, b::Int)
#     s = it.siteinds(psi)  
#     it.orthogonalize!(psi, b)
#     _,S = it.svd(psi[b], (it.linkind(psi, b-1), s[b]))
#     SvN = 0.0
#     for n in 1:it.dim(S, 1)
#       p = S[n,n]^2
#       SvN -= p * log(p)
#     end
#     return SvN
# end

function entanglement_entropy(psi::it.MPS) #fix
    #pastaq
    ψ = la.normalize!(copy(psi))
    N = length(ψ)
    bond = N ÷ 2
    it.orthogonalize!(ψ, bond)
  
    row_inds = (it.linkind(ψ, bond - 1), it.siteind(ψ, bond))
    u, s, v = it.svd(ψ[bond], row_inds)
  
    S = 0.0
    for n in 1:it.dim(s, 1)
      λ = s[n, n]^2
      S -= λ * log(λ + 1e-20)
    end
    return S
end