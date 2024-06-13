"""
_get_codeword
"""
get_codeword(codestate::sa.SparseVector)=fock_basis.(findall(x->!isapprox(x,0;atol=1000eps()),codestate).-1,get_N(codestate))


function stabilizers_to_generator(stabilizers::Vector)

    num_qubits = length(split(stabilizers[1],",")) #this is n
    num_stabilizers = length(stabilizers)
    
    generator_matrix = zeros(Int, num_stabilizers, 2 * num_qubits)
    
    for (i, stabilizer) in enumerate(stabilizers)
        for (j, op) in enumerate(split(stabilizer,","))
            if op == "X"
                generator_matrix[i, j] = 1
            elseif op == "Z"
                generator_matrix[i, j + num_qubits] = 1
            elseif op == "Y"
                generator_matrix[i, j] = 1  # X part
                generator_matrix[i, j + num_qubits] = 1  # Z part
            end
        end
    end
    
    return generator_matrix
end

"""
standard_form_logicals(G::AbstractMatrix)
"""
function standard_form_logicals(G::AbstractMatrix)
    n, m = size(G, 2) ÷ 2, size(G, 1)
    k = n - m
    
    generator = G[:, 1:n]
    G2 = G[:, (n+1):end]
    
    # Step 1: Gaussian Elimination on generator & adjust G2 rows & cols accordingly
    generator_rref, r, generator_transform_rows, generator_transform_cols = reduced_row_echelon(generator)
    G = hcat(generator_rref, mod.(generator_transform_rows * G2 * generator_transform_cols, 2))
    
    # Step 2: Gaussian Elimination on E
    A = G[1:r, (r+1):n]
    B = G[1:r, n+1:(n+r)]
    C = G[1:r, (n+r+1):end]
    D = G[(r+1):end, n+1:(n+r)]
    E = G[(r+1):end, (n+r+1):end]
    
    if n - k - r != 0
        E_rref, s, E_transform_rows, E_transform_cols = reduced_row_echelon(E)
        
        A = mod.(A * E_transform_cols, 2)
        C = mod.(C * E_transform_cols, 2)
        D = mod.(E_transform_rows * D, 2)
        
        G_row_1 = hcat(hcat(Matrix{Int}(la.I, r, r), A), hcat(B, C))
        G_row_2 = hcat(zeros(Int, n-k-r, n), hcat(D, E_rref))
        G_standard = vcat(G_row_1, G_row_2)
    end

    id, A1, A2, B, C1, C2, D, E2 = _extract_blocks(G_standard,r)

    logical_Z_vec=[fill(0,n)...,A2'...,fill(0,m-r)...,fill(1,k)...]
    logical_X_vec=[fill(0,r)...,E2'...,fill(1,k)...,E2'C1'+C2'...,fill(0,m-r)...,fill(0,k)...]

    return G_standard, logical_vec_2_ops.([logical_X_vec,logical_Z_vec]), (logical_X_vec,logical_Z_vec)
end

function logical_vec_2_ops(logical_vec::Vector)

    n = length(logical_vec) ÷ 2
    
    ops=[]
    for (i,x)=enumerate(logical_vec)
        if x==1
            if i<=n
                o=logical_vec[i+n] == 1 ? Op("Y",i) : Op("X",i)
                push!(ops,o)
            elseif i>n
                o=logical_vec[i-n] == 1 ? Op("Y",i-n) : Op("Z",i-n)
                push!(ops,o)
            end
        end
    end
    
    return ops
    end



function _extract_blocks(G_standard_form,r::Int)
    #block sizes: r, m-r, k
    # r=la.rank(G);println("rank of G=$(r)")
    n, m = size(G_standard_form, 2) ÷ 2, size(G_standard_form, 1)
    k = n - m

    left=r+1
    block_size=(m-r)
    A1 = G_standard_form[1:r, left:left+block_size-1]

    left+=block_size
    block_size=k
    A2 = G_standard_form[1:r, left:left+block_size-1]

    left+=block_size
    block_size=r
    B = G_standard_form[1:r, left:left+block_size-1]

    left+=block_size
    block_size=(m-r)
    C1 = G_standard_form[1:r, left:left+block_size-1]

    left+=block_size
    block_size=k
    C2 = G_standard_form[1:r, left:left+block_size-1]

    #row2

    left=n+1
    block_size=r
    D = G_standard_form[(r+1):end, left:left+block_size-1]

    left+=block_size
    block_size=(m-r)
    id = G_standard_form[(r+1):end, left:left+block_size-1]

    left+=block_size
    block_size=k
    E2 = G_standard_form[(r+1):end, left:left+block_size-1]

    return id, A1, A2, B, C1, C2, D, E2
end


"""
    StabilizerCode(n::Int,k::Int,d::Int,stabilizers::Vector,logicals::Dict)

    n=4
    k=2
    d=2

    logicals=Dict()

    logicals["I",1]=[("I",1)]
    logicals["I",2]=[("I",1)]

    logicals["X",1]=[("X",1),("X",2)]
    logicals["X",2]=[("X",1),("X",4)]
    logicals["Z",1]=[("Z",1),("Z",4)]
    logicals["Z",2]=[("Z",1),("Z",2)]
    logicals["CNOT",1,2]=[("SWAP",3,4)]
    logicals["CNOT",2,1]=[("SWAP",3,4),("SWAP",2,4),("SWAP",1,2)]

    stabilizers=[
        ("X,X,X,X"),
        ("Z,Z,Z,Z")
    ]

    code422=StabilizerCode(n,k,d,stabilizers,logicals)

    fields(code422)
"""
struct StabilizerCode #alpha version
    n::Int
    k::Int
    d::Int
    stabilizers::Vector
    generator::AbstractMatrix
    generator_standard::AbstractMatrix
    logical_XZ_ops::Vector
    logicals::Dict
    codestates::Vector
    codewords::Vector
    stabilizer_matrix::Vector
    stats::NamedTuple
    ops::Function
    apply::Function
    encode::Function
    decode::Function
    syndrome::Function
    correct::Function
    circuit::Function

    function StabilizerCode(n::Int,k::Int,d::Int,stabilizers::Vector{String},logicals::Dict)

        # if length(first(logicals)) != k
            # throw("k and logicals are incompatible")
        if length(split(stabilizers[1],","))!=n
            throw("n and stabilizers are incompatible")
        end

        ##find init from stabilizers
        stabilizer_matrix=string_to_matrix.(stabilizers)
        generator=stabilizers_to_generator(stabilizers)
        generator_standard, logical_XZ_ops, XZ_vecs = standard_form_logicals(generator)
        len_stabilizers=length(stabilizers)
        e,v=la.eigen(Matrix(sum(stabilizer_matrix)/len_stabilizers))
        codestates = [sa.sparse(v[:, i]) for i in 1:length(e) if abs(e[i] - 1) < 1000eps()];sa.droptol!.(codestates,1000eps());
        len_codestates=length(codestates)
        codewords=[get_codeword(code_s) for code_s=codestates]
        if length(codestates) != 2^k
            throw("wrong stabilizers")
        end

        ##

        function new_encode(state::sa.SparseVector)

            throw("fix")
        end

        function new_decode(state::sa.SparseVector)

            throw("fix")
        end

        function new_syndrome(state::sa.SparseVector)

            throw("fix")
        end

        function new_correct(state::sa.SparseVector)

            throw("fix")
        end

        function new_circuit(state::sa.SparseVector)

            throw("fix")
        end

        new_ops(name::String, logical_qubit::Int)=[Op(x...) for x=logicals[name,logical_qubit]]
        new_ops(name::String, logical_qubit::Int, target_qubit::Int)=[Op(x...) for x=logicals[name,logical_qubit,target_qubit]]

        function new_ops(op::Op)
            if op.q == 1
                return [Op(x...) for x=logicals[op.name,op.qubit]]
            elseif op.q == 2
                return [Op(x...) for x=logicals[op.name,op.qubit,op.target_qubit]]
            end
        end

        #one qubit
        function new_apply(state::Union{sa.SparseVector,sa.SparseMatrixCSC}, name::String, qubit::Int;noise::Union{Bool,NoiseModel}=false)
            ops=new_ops(name, qubit)
            N=get_N(state)

            if qubit > N
                throw("qubit cannot be larger than N=$(N)")
            end

            for o=ops
                state=apply(state,o;noise=noise)
            end

            return state
        end

        #two qubit
        function new_apply(state::Union{sa.SparseVector,sa.SparseMatrixCSC}, name::String, logical_qubit::Int, logical_target::Int;noise::Union{Bool,NoiseModel}=false)
            N=get_N(state)
            max_qt=max(logical_qubit,logical_target)
            min_qt=min(logical_qubit,logical_target)
            distance=abs(logical_qubit-logical_target)

            if max_qt > N
                throw("qubit or target_qubit cannot be larger than N=$(N)")
            elseif distance==0
                throw("qubit and target cannot be same")
            elseif abs(logical_qubit-logical_target)>1
                throw("nonlocal logical operations are not supported")
            end

            keyexist=(name, logical_qubit, logical_target) ∈ keys(logicals)

            if keyexist
                ops=new_ops(name, logical_qubit, logical_target)
            else
                throw("logical does not exist")
            end

            for o=ops
                state=apply(state,o;noise=noise)
            end

            return state
        end

        new_apply(state::Union{sa.SparseVector,sa.SparseMatrixCSC}, op::Op; noise::Union{Bool,NoiseModel}=false)=op.q==1 ? new_apply(state, op.name, op.qubit;noise=noise) : new_apply(state, op.name, op.qubit, op.target_qubit;noise=noise)
        
        stats=(
            len_logical=length(keys(logicals)),
            len_stabilizers=len_stabilizers,
            len_codestates=len_codestates,
            XZ_vecs=XZ_vecs
        )

        return new(n,k,d,stabilizers,generator,generator_standard,logical_XZ_ops,logicals,codestates,codewords,stabilizer_matrix,stats,new_ops,new_apply,new_encode,new_decode,new_syndrome,new_correct,new_circuit);

    end

end