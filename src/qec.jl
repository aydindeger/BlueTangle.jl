"""
get_codeword(codestate::sa.SparseVector)
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
get_standard_form(G::AbstractMatrix)
"""
function get_standard_form(G::AbstractMatrix)
    n, m = size(G, 2) ÷ 2, size(G, 1)
    k = n - m
    
    generator = G[:, 1:n]
    G2 = G[:, (n+1):end]

    generator_rref, r, generator_transform_rows, generator_transform_cols = reduced_row_echelon(generator)
    G = hcat(generator_rref, mod.(generator_transform_rows * G2 * generator_transform_cols, 2))
    
    # Gaussian Elimination on E
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

    return G_standard, r
end


function get_XZ_logicals(G_standard::AbstractMatrix,logical_XZ_ops::Dict,r::Int)

    n, m = size(G_standard, 2) ÷ 2, size(G_standard, 1)
    k = n - m

    id, A1, A2, B, C1, C2, D, E = _extract_blocks(G_standard,r)

    logical_XZ_vecs=Dict()

    if k>1
        println("something is wrong with logicals for k>1")
    end

    U0=zeros(k,r)
    U2=E'
    U3=la.I(k)
    V1=(E'*C1')+C2'
    U1=zeros(k,m-r)
    V3=zeros(k,k)
    #size: r m-r k
    X_vecs=[U0 U2 U3 V1 U1 V3]
    Z_vecs=[zeros(k,n) A2' U1 la.I(k)]

    for ki=1:k
        logical_XZ_vecs["logZ",ki]=Z_vecs[ki,:]
        logical_XZ_vecs["logX",ki]=X_vecs[ki,:]
        logical_XZ_ops["logZ",ki]=logical_vec_2_ops(Z_vecs[ki,:])
        logical_XZ_ops["logX",ki]=logical_vec_2_ops(X_vecs[ki,:])
    end

    return logical_XZ_ops, logical_XZ_vecs

end


function logical_vec_2_ops(logical_vec::Vector)

    n = length(logical_vec) ÷ 2
    
    ops=[]
    for (i,x)=enumerate(logical_vec)
        if x==1
            if i<=n
                o=logical_vec[i+n] == 1 ? ("Y",i) : ("X",i)
                push!(ops,o)
            elseif i>n
                o=logical_vec[i-n] == 1 ? ("Y",i-n) : ("Z",i-n)
                push!(ops,o)
            end
        end
    end
    
    return ops
    end

function _extract_blocks(G_standard_form,r::Int) #r=rank
    #block sizes: r, m-r, k
    n, m = size(G_standard_form, 2) ÷ 2, size(G_standard_form, 1)
    k = n - m

    id=G_standard_form[1:r,1:r]
    A1=G_standard_form[1:r,r+1:n-k]
    A2=G_standard_form[1:r,n-k+1:n]
    B=G_standard_form[1:r,n+1:n+r]
    C1=G_standard_form[1:r,n+r+1:2n-k]
    C2=G_standard_form[1:r,2n-k+1:2n]
    
    D=G_standard_form[r+1:m,n+1:n+r]
    id2=G_standard_form[r+1:m,n+r+1:2n-k]
    E=G_standard_form[r+1:m,2n-k+1:2n]

    return id, A1, A2, B, C1, C2, D, E
end

function encoding_circuit_from_generator(generator_standard::AbstractMatrix,logical_XZ_vecs::Dict,r::Int)

    #stac ebook, gottesman's thesis
    m,n=size(generator_standard,1),size(generator_standard,2) ÷ 2
    k=n-m

    standard_generators_x=generator_standard[:,1:n]
    standard_generators_z=generator_standard[:,n+1:end]

    encoding_circuit = Vector{Op}()

    # First loop: Add CX gates based on logical_xs matrix
    for i in 1:k
        for j in r+1:n-k
            if logical_XZ_vecs["logX",i][j]==1
                push!(encoding_circuit, Op("CX", n-k+i, j))
            end
        end
    end

    # Second loop: Add H, CX, and CZ gates based on standard_generators_x and standard_generators_z matrices
    for i in 1:r
        push!(encoding_circuit, Op("H", i))
        for j in 1:n
            if i == j
                continue
            end
            if (standard_generators_x[i, j]==1) && (standard_generators_z[i, j]==1)
                push!(encoding_circuit, Op("CX", i, j))
                push!(encoding_circuit, Op("CZ", i, j))
            elseif standard_generators_x[i, j]==1
                push!(encoding_circuit, Op("CX", i, j))
            elseif standard_generators_z[i, j]==1
                push!(encoding_circuit, Op("CZ", i, j))
            end
        end
    end

    return encoding_circuit
end


"""
    StabilizerCode(stabilizers::Vector,logicals::Dict)

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

    code422=StabilizerCode(stabilizers,logicals)

    fields(code422)
"""
struct StabilizerCode #alpha version
    n::Int
    k::Int
    d::Int
    stabilizers::Vector
    generator::AbstractMatrix
    generator_standard::AbstractMatrix
    logicals::Dict
    codestates::Vector
    codewords::Vector
    encoding_ops::Vector
    stabilizer_matrix::Vector
    info::NamedTuple
    ops::Function
    apply::Function
    encode::Function
    decode::Function
    syndrome::Function
    correct::Function
    circuit::Function

    function StabilizerCode(stabilizers::Vector{String},logicals::Dict;d::Int=-1)

        n=length(split.(stabilizers,",")[1])
        m=length(stabilizers)
        k=n-m

        ##find init from stabilizers
        stabilizer_matrix=string_to_matrix.(stabilizers)
        generator=stabilizers_to_generator(stabilizers)
        generator_standard, r = get_standard_form(generator)
        logicals, logical_XZ_vecs = get_XZ_logicals(generator_standard,logicals,r)
        encoding_ops=encoding_circuit_from_generator(generator_standard,logical_XZ_vecs,r)
        e,v=la.eigen(Matrix(sum(stabilizer_matrix)/m))
        codestates = [-sa.sparse(v[:, i]) for i in 1:length(e) if abs(e[i] - 1) < 1000eps()];sa.droptol!.(codestates,1000eps());
        len_codestates=length(codestates)
        codewords=[get_codeword(code_s) for code_s=codestates]
        if length(codestates) != 2^k
            throw("wrong stabilizers")
        end

        ##

        function new_encode(state_init::sa.SparseVector)

            if k>1
                println("something is wrong with logicals for k>1")
            end        

            N=get_N(state_init)
            if N!=k
                throw(ArgumentError("init state cannot be larger than k=$(k) of the code."))
            end

            state=kron(zero_state(n-k),state_init)

            for o=encoding_ops
                state=apply(state,o)
            end

            return state

        end

        function new_decode(state::sa.SparseVector)

            for o=encoding_ops'
                state=apply(state,o)
            end

            return sa.sparse(sa.findnz(state)[2])

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
        
        info=(
            len_logical=length(keys(logicals)),
            len_codestates=len_codestates,
            logical_XZ_vecs=logical_XZ_vecs,
            r=r,
            m=m #len_stabilizers
        )

        return new(n,k,d,stabilizers,generator,generator_standard,logicals,codestates,codewords,encoding_ops,stabilizer_matrix,info,new_ops,new_apply,new_encode,new_decode,new_syndrome,new_correct,new_circuit);

    end

end