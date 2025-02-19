"""
get_codeword(codestate::AbstractVectorS)
"""
get_codeword(codestate::AbstractVectorS)=fock_basis.(findall(x->!isapprox(x,0;atol=1000eps()),codestate).-1,get_N(codestate))

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
    
    E_rref, s, E_transform_rows, E_transform_cols = reduced_row_echelon(E)

    if n - k - r != 0
        A = mod.(A * E_transform_cols, 2)
        C = mod.(C * E_transform_cols, 2)
        D = mod.(E_transform_rows * D, 2)
        
        G_row_1 = hcat(hcat(Matrix{Int}(la.I, r, r), A), hcat(B, C))
        G_row_2 = hcat(zeros(Int, n-k-r, n), hcat(D, E_rref))
        G = vcat(G_row_1, G_row_2)
    end
    
    permute_matrix_E = Matrix{Int64}(la.I, n, n)
    len_E=size(E_transform_cols,1)
    permute_matrix_E[n-len_E+1:n, n-len_E+1:n] .= E_transform_cols
    permute_matrix_inv=permute_matrix_E'*generator_transform_cols'

    if r!=rank_of_rref(G)
        throw("something's wrong with the rank or standard form")
    end

    permute_vector_inv = [findfirst(permute_matrix_inv'[:, i] .== 1) for i in 1:n] #' is correct
    logical_indices=[permute_vector_inv[n-k+i] for i=1:k]

    return G, permute_matrix_inv, permute_vector_inv, logical_indices
end

"""
get_XZ_logicals(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix)
"""
function get_XZ_logicals(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix)
    get_XZ_logicals!(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix,Dict())
end

"""
get_XZ_logicals!(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix,logical_XZ_ops::Dict)
"""
function get_XZ_logicals!(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix,logical_XZ_ops::Dict)

    r=rank_of_rref(G_standard)
    n, m = size(G_standard, 2) ÷ 2, size(G_standard, 1)
    k = n - m

    id, A1, A2, B, C1, C2, D, E = _extract_blocks(G_standard,r)

    logical_XZ_vecs=Dict()

    U0=zeros(k,r)
    U2=E'
    U3=la.I(k)
    V1=(E'*C1')+C2'
    U1=zeros(k,m-r)
    V3=zeros(k,k)
    #size: r m-r k
    X_vecs=[U0 U2 U3 V1 U1 V3]
    Z_vecs=[zeros(k,n) A2' U1 la.I(k)]

    # #permute
    ip=permute_matrix_inv#reduced_row_echelon_inverse(permute_matrix)
    Xvecs=hcat(X_vecs[:,1:n]*ip,X_vecs[:,n+1:end]*ip)
    Zvecs=hcat(Z_vecs[:,1:n]*ip,Z_vecs[:,n+1:end]*ip)

    for ki=1:k
        logical_XZ_vecs["Z",ki]=sa.sparse(Zvecs[ki,:])
        logical_XZ_vecs["X",ki]=sa.sparse(Xvecs[ki,:])
        logical_XZ_ops["Z",ki]=logical_vec_2_ops(Zvecs[ki,:])
        logical_XZ_ops["X",ki]=logical_vec_2_ops(Xvecs[ki,:])
    end

    return logical_XZ_ops, sa.sparse(Xvecs), sa.sparse(Zvecs)

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

    if (la.one(id) ≈ id) || (la.one(id2) ≈ id2)
        return id, A1, A2, B, C1, C2, D, E
    else
        throw("something's wrong with the standard form")
    end

end

"""
encoding_circuit_from_generator(generator_standard::AbstractMatrix,Xvecs::AbstractMatrixS,permute_vector_inv::AbstractVector)
"""
function encoding_circuit_from_generator(generator_standard::AbstractMatrix,Xvecs::AbstractMatrixS,permute_vector_inv::AbstractVector)

    m,n=size(generator_standard,1),size(generator_standard,2) ÷ 2
    r=rank_of_rref(generator_standard)
    k=n-m
    
    encoding_circuit = Vector{Op}()
    
    # First loop: Add CX gates based on logical_xs matrix
    for i in 1:k #goes over whole
    
        control_qubit=permute_vector_inv[n-k+i]
    
        # if Xvecs[i,control_qubit]==1
        #    #control qubit 
        # end
    
        for j in 1:n
            if control_qubit != j
                if Xvecs[i,j]==1
                    push!(encoding_circuit,Op("CX",control_qubit,j))
                end
            end
        end
    
    end
    
    for i in 1:r
    
        control_qubit=i
        if generator_standard[i,i+n]==0
            push!(encoding_circuit, Op("H", i))
        else
            push!(encoding_circuit, Op("H", i))
            push!(encoding_circuit, Op("S", i))
        end
    
        for j in 1:n
            if i != j
                if (generator_standard[i, j]==1) && (generator_standard[i, j+n]==0)
                    push!(encoding_circuit, Op("CX", control_qubit,j))
                elseif (generator_standard[i, j]==0) && (generator_standard[i, j+n]==1)
                    push!(encoding_circuit, Op("CZ", control_qubit, j))
                elseif (generator_standard[i, j]==1) && (generator_standard[i, j+n]==1)
                    push!(encoding_circuit, Op("CY", control_qubit,j))
                end
            end
        end
    end
    
    return encoding_circuit
end

function get_ops_syndrome(generator_standard::AbstractMatrix)

    m,n=size(generator_standard,1),size(generator_standard,2) ÷ 2

    standard_generators_x=generator_standard[:,1:n]
    standard_generators_z=generator_standard[:,n+1:end]

    syndrome_circuit = Vector{Op}()

    for i in 1:m
        push!(syndrome_circuit, Op("H", n + i))
    end

    for i in 1:m
        for j in 1:n
            if (standard_generators_x[i, j]==1) && (standard_generators_z[i, j]==1)
                push!(syndrome_circuit, Op("CX", n + i, j))
                push!(syndrome_circuit, Op("CZ", n + i, j))
            elseif standard_generators_x[i, j]==1
                push!(syndrome_circuit, Op("CX", n + i, j))
            elseif standard_generators_z[i, j]==1
                push!(syndrome_circuit, Op("CZ", n + i, j))
            end
        end
    end

    for i in 1:m
        push!(syndrome_circuit, Op("H", n + i))
    end

    for i in 1:m
        push!(syndrome_circuit, Op("MZ", n + i))
    end

    return syndrome_circuit
end

"""
    identify_error(syndrome::AbstractVector, generator_standard)

    Simple decoder will return the erroneous op
"""
function identify_error(syndrome::AbstractVector, generator_standard::AbstractMatrix)

    if length(syndrome) != size(generator_standard, 1)
        throw("Syndrome length mismatch: expected $(size(generator_standard, 1)), got $(length(syndrome))")
    end

    n = size(generator_standard, 2) ÷ 2
    # Iterate through all possible single qubit errors
    for qubit in 1:n
        for error_op in ["X", "Z", "Y"]
            # Create an error vector
            error_vec = zeros(Int, 2 * n)
            error_vec[qubit] = (error_op in ["Z", "Y"]) * 1  # Set Z or Y error
            error_vec[n + qubit] = (error_op in ["X", "Y"]) * 1  # Set X or Y error
            
            # Calculate the syndrome for this error_op
            calculated_syndrome = mod.(generator_standard * error_vec, 2)
            
            # Check if the calculated syndrome matches the provided syndrome
            if all(calculated_syndrome .== syndrome)
                println("Error $(error_op) on qubit $(qubit)")
                return Op(error_op,qubit)
            end
        end
    end
    return "No single-qubit error matches this syndrome"
end

"""
physical_ops(name::String, logical_qubit::Int, logicals::Dict)
"""
physical_ops(name::String, logical_qubit::Int, logicals::Dict)=[Op(x...) for x=logicals[name,logical_qubit]]

"""
physical_ops(name::String, logical_qubit::Int, target_qubit::Int, logicals::Dict)
"""
physical_ops(name::String, logical_qubit::Int, target_qubit::Int, logicals::Dict)=[Op(x...) for x=logicals[name,logical_qubit,target_qubit]]

"""
physical_ops(op::QuantumOps, logicals::Dict)
"""
function physical_ops(op::QuantumOps, logicals::Dict)
    if isa(op,OpF)
        vals=logicals[op]
    elseif op.q == 1
        vals=logicals[op.name,op.qubit]
    elseif op.q == 2
        vals=logicals[op.name,op.qubit,op.target_qubit]
    end

    return isa(vals,QuantumOps) ? vals : (isa(first(vals),QuantumOps) ? vals : [Op(x...) for x=vals])
end

"""
physical_ops(ops::Vector{QuantumOps}, logicals::Dict)
"""
function physical_ops(ops::Vector{QuantumOps}, logicals::Dict)
    list_of_ops=Vector{QuantumOps}()
    for op=ops
        append!(list_of_ops,physical_ops(op, logicals))
    end

    return list_of_ops
end

function logical_string_to_dict(logical_string::String)
    lines = split(chomp(logical_string), '\n')
    dict = Dict()
    for (line_index, line) in enumerate(lines)
        for (char_index, char) in enumerate(collect(line))
            if char != 'I'
                key = (char, line_index)
                values = [(char, pos) for (pos, c) in enumerate(collect(line)) if c == char]
                dict[key] = values
            end
        end
    end
    return dict
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

    stabilizers=[
        ("X,X,X,X,I,I,I"),
        ("X,X,I,I,X,X,I"),
        ("X,I,X,I,X,I,X"),
        ("Z,Z,Z,Z,I,I,I"),
        ("Z,Z,I,I,Z,Z,I"),
        ("Z,I,Z,I,Z,I,Z")
    ]

    code=StabilizerCode(stabilizers,logicals)

    fields(code)


    EXAMPLES:

    #this encodes and decodes

    a=random_state(code.k)
    b=code.decode(code.encode(a))
    #a and b are same up to a global phase
    expect(a,"Z")[1] ≈ expect(b,"Z")[1]

    #error correction

    init_state=zero_state(code.k)
    code_zero=code.encode(init_state)
    code_zero_err=apply(code_zero,Op("Y",1)) #error on first qubit
    code_zero_ancilla,syndrome,op=code.syndrome(code_zero_err)
    code_final=code.correct(code_zero_ancilla,syndrome)
    code_decoded=code.decode(code_final)

    expect(init_state,"Z")[1] ≈ expect(code_decoded,"Z")[1]
"""
struct StabilizerCode #alpha version
    n::Int
    k::Int
    d::Int
    m::Int
    stabilizers::Vector
    generator::AbstractMatrix
    generator_standard::AbstractMatrix
    permute_matrix_inv::AbstractMatrix
    permute_vector_inv::AbstractVectorS
    logical_indices::AbstractVector
    logicals::Dict
    # codestates::Vector
    # codewords::Vector
    ops_encoding::Vector #encoding circuit
    ops_syndrome::Vector #syndrome circuit
    # stabilizer_matrix::Vector
    info::NamedTuple
    ops::Function
    apply::Function
    encode::Function
    decode::Function
    syndrome::Function
    correct::Function

    function StabilizerCode(stabilizers::Union{Vector{String},String},logicals::Dict;d::Int=-1)

        if isa(stabilizers,String)
            stabilizers = [join(collect(line), ",") for line in split(chomp(stabilizers), '\n')]
        end

        n=length(split.(stabilizers,",")[1])
        m=length(stabilizers)
        k=n-m

        ##find init from stabilizers
        # stabilizer_matrix=string_to_matrix.(stabilizers)
        generator=stabilizers_to_generator(stabilizers)
 
        generator_standard, permute_matrix_inv, permute_vector_inv, logical_indices = get_standard_form(generator)
        logicals, Xvecs, Zvecs = get_XZ_logicals!(generator_standard,permute_matrix_inv,logicals)
        ops_encoding=encoding_circuit_from_generator(generator_standard,Xvecs,permute_vector_inv)
        ops_syndrome=get_ops_syndrome(generator_standard)

        # below for codespace and codewords
        # e,v=la.eigen(Matrix(sum(stabilizer_matrix)/m))
        # codestates = [sa.sparse(v[:, i]) for i in 1:length(e) if abs(e[i] - 1) < 1000eps()];sa.droptol!.(codestates,1000eps());
        # len_codestates=length(codestates)
        # codewords=[get_codeword(code_s) for code_s=codestates]
        # if length(codestates) != 2^k
        #     throw("wrong stabilizers")
        # end

        ##

        new_ops(name::String, logical_qubit::Int)=physical_ops(name,logical_qubit,logicals)
        new_ops(name::String, logical_qubit::Int, target_qubit::Int)=physical_ops(name,logical_qubit,target_qubit,logicals)
        new_ops(op::QuantumOps)=physical_ops(op,logicals)
        new_ops(ops::Vector{QuantumOps})=physical_ops(ops,logicals)

        # function new_encode(state_init::AbstractVectorS;noise::Union{Bool,NoiseModel}=false,encoding::Vector=[]) #encoded state is last  

        # new_encode= :plus, :minus, :T or :zero
        
        new_encode(state_init_sym::Symbol=:zero;random_op_count::Int=20,noise::Union{Bool,NoiseModel}=false,encoding::Vector=[])=new_encode(n,state_init_sym;random_op_count=random_op_count,noise=noise,encoding=encoding)

        #mps
        function new_encode(N::Union{Int,Vector{it.Index{Int64}}},state_init_sym::Symbol=:zero;random_op_count::Int=20,noise::Union{Bool,NoiseModel}=false,encoding::Vector=[])

            state=qec_state_prep(N,logical_indices,state_init_sym;random_op_count=random_op_count)[1]#[1] indicates qec state

            ops=isempty(encoding) ? ops_encoding : encoding
            for o=ops
                state=apply(state,o;noise=noise)
            end

            return state

        end

        function new_decode(state::AbstractVectorS;noise::Union{Bool,NoiseModel}=false,decoding::Vector=[])

            ops=isempty(decoding) ? ops_encoding' : decoding

            for o=ops
                state=apply(state,o;noise=noise)
            end

            return reduce_and_relabel_qubits(state, logical_indices)

        end

        function new_syndrome(state_encoded::AbstractVectorS;noise::Union{Bool,NoiseModel}=false)

            throw("fix: needs to be fixed")

            if get_N(state_encoded)!=n
                throw("input state must be encoded in the code")
            end

            state_encoded_ancillas=la.kron(state_encoded,zero_state(m)) #add m ancillas to the end

            for o=ops_syndrome
                state_encoded_ancillas=apply(state_encoded_ancillas,o;noise=noise)
            end

            syndrome_vec=sample_bit(state_encoded_ancillas,1)[1][n+1:end]
            # syndrome_pos=findall(!=(0),syndrome_vec)

            if sum(syndrome_vec)==0 #isempty(syndrome_pos)
                println("No errors detected!")
                op_error=Op("I",1;noisy=false)
            else
                op_error=identify_error(syndrome_vec, generator_standard)
            end

            for i=n+1:n+m #reset ancilla
                state_encoded_ancillas=Op("RES",i)*state_encoded_ancillas
            end

            return state_encoded_ancillas,syndrome_vec,op_error
        end

        function new_correct(state_encoded_ancillas::AbstractVectorS,syndrome_vec::Vector;noise::Union{Bool,NoiseModel}=false)

            throw("fix: needs to be fixed")
            if 0<d<=2
                throw("The code cannot correct errors!")
            end
                
            if sum(syndrome_vec)==0 #isempty(syndrome_pos)
                println("No errors detected!")
                return state_encoded_ancillas
            end

            op_error=identify_error(syndrome_vec, generator_standard)
            final_state=apply(state_encoded_ancillas,op_error;noise=noise)
            print("error is corrected!")

            return final_state

        end

        function new_correct(state_encoded::AbstractVectorS;noise::Union{Bool,NoiseModel}=false) #firstly measure syndromes

            throw("fix: needs to be fixed")
            if 0<d<=2
                throw("The code cannot correct errors!")
            end
            
            state_encoded_ancillas,syndrome_vec,op_error=new_syndrome(state_encoded::AbstractVectorS)
            
            if sum(syndrome_vec)==0 #isempty(syndrome_pos)
                println("No errors detected!")
                return state_encoded_ancillas
            end

            final_state=apply(state_encoded_ancillas,op_error;noise=noise)
            print("error is corrected!")

            return final_state

        end

        #one qubit
        function new_apply(state::AbstractVectorS, name::String, qubit::Int;noise::Union{Bool,NoiseModel}=false)
            ops=physical_ops(name, qubit, logicals)
            N=get_N(state)

            if qubit > N+m
                throw("qubit cannot be larger than system+ancilla qubits: $(N+m)")
            end

            for o=ops
                state=apply(state,o;noise=noise)
            end

            return state
        end

        #two qubit
        function new_apply(state::AbstractVectorS, name::String, logical_qubit::Int, target_qubit::Int;noise::Union{Bool,NoiseModel}=false)
            N=get_N(state)
            max_qt=max(logical_qubit,target_qubit)
            min_qt=min(logical_qubit,target_qubit)
            distance=abs(logical_qubit-target_qubit)

            if max_qt > N+m
                throw("qubit or target_qubit cannot be larger than system+ancilla qubits: $(N+m)")
            elseif distance==0
                throw("qubit and target cannot be same")
            elseif abs(logical_qubit-target_qubit)>1
                throw("nonlocal logical operations are not supported")
            end

            keyexist=(name, logical_qubit, target_qubit) ∈ keys(logicals)

            if keyexist
                ops=physical_ops(name,logical_qubit,target_qubit,logicals)
            else
                throw("logical does not exist")
            end

            for o=ops
                state=apply(state,o;noise=noise)
            end

            return state
        end

        new_apply(state::AbstractVectorS, op::Op; noise::Union{Bool,NoiseModel}=false)=op.q==1 ? new_apply(state, op.name, op.qubit;noise=noise) : new_apply(state, op.name, op.qubit, op.target_qubit;noise=noise)
        
        info=(
            permute_vector=[findfirst(permute_matrix_inv[:, i] .== 1) for i in 1:n],
            len_logical=length(keys(logicals)),
            # len_codestates=len_codestates,
            logical_XZ_vecs=(Xvecs,Zvecs)
        )

        return new(n,k,d,m,stabilizers,generator,generator_standard,permute_matrix_inv,permute_vector_inv,logical_indices,logicals,ops_encoding,ops_syndrome,info,new_ops,new_apply,new_encode,new_decode,new_syndrome,new_correct);

    end

end

function permutation_vector_from_matrix(permute_matrix::AbstractMatrix{Int})
    n = size(permute_matrix, 1)
    p = Vector{Int}(undef, n)
    for i in 1:n
        p[i] = findfirst(permute_matrix[:, i] .== 1)
    end
    return p
end


"""
qec_state_prep(n::Union{Int,Vector{it.Index{Int64}}},logical_indices::Vector,state_init_sym::Symbol=:zero;random_op_count::Int=20)
"""
function qec_state_prep(n::Union{Int,Vector{it.Index{Int64}}},logical_indices::Vector,state_init_sym::Union{Symbol,Vector}=:zero;random_op_count::Int=20)
    # state_init_sym=:random

    mps_bool=typeof(n)==Int ? false : true
    
    len_k=length(logical_indices)
    k=mps_bool==false ? len_k : n[1:len_k]

    state=zero_state(n)
    state_logical=zero_state(k)

    if (state_init_sym==:zeros) || (state_init_sym==:zero)
        return state,state_logical
    elseif (state_init_sym==:ones) || (state_init_sym==:one)
        state=one_state(n)
        state_logical=one_state(k)
    elseif state_init_sym==:plus
        for i=logical_indices
            state=Op("H",i)*state
        end
        for i=1:len_k
            state_logical=Op("H",i)*state_logical
        end
    elseif state_init_sym==:minus
        state=one_state(n)
        state_logical=one_state(k)
        for i=logical_indices
            state=Op("H",i)*state
        end
        for i=1:len_k
            state_logical=Op("H",i)*state_logical
        end
    elseif state_init_sym==:random
        ops_random=BlueTangle.random_ops(len_k,random_op_count)
        for o=ops_random     
            state_logical=o*state_logical
            if o.q==1
                state=Op(o.name,logical_indices[o.qubit])*state
            elseif o.q==2
                state=Op(o.name,logical_indices[o.qubit],logical_indices[o.target_qubit])*state
            end
        end
    elseif isa(state_init_sym,AbstractVector)
        state_logical=mps_bool==true ? product_state(k,state_init_sym) : product_state(state_init_sym)

        state_fock=mps_bool==true ? Int.(la.zeros(length(n))) : Int.(la.zeros(n))
        state_fock[logical_indices] .= state_init_sym

        state=mps_bool==true ? product_state(n,state_fock) : product_state(state_fock)
    else
        throw("error: invalid state name!")
    end

    return state,state_logical
end

"""
reduce_and_relabel_qubits(psi::AbstractVectorS, logical_indices::Vector,qubit_mapping::Vector=Vector[])
logical_indices=qubit_mapping
"""
function reduce_and_relabel_qubits(psi::AbstractVectorS, logical_indices::Vector,qubit_mapping::Vector=Vector[])
    n = Int(round(log2(length(psi))))  # Total number of qubits
    k = length(logical_indices)          # Number of qubits to keep
    psi_new = sa.spzeros(ComplexF64, 2^k)   # New state vector for k qubits

    if qubit_mapping!=Vector[]
        logical_indices = [findfirst(x->x==y,qubit_mapping) for y=logical_indices]
    end
    
    # Iterate over all basis states in the reduced Hilbert space
    for j in 0:(2^k - 1)
        # Binary representation of the new basis state (reversed due to endianness)
        b_new = reverse(digits(j, base=2, pad=k))
        
        # Build the corresponding full basis state in the original Hilbert space
        b_full = zeros(Int, n)
        for idx in 1:k
            b_full[logical_indices[idx]] = b_new[idx]
        end
        # The discarded qubits are assumed to be in |0⟩ state (bits set to 0)
        
        # Convert the full binary string to an index in the original state vector
        i = 0
        for bit in b_full
            i = (i << 1) | bit  # Left shift and add the bit
        end
        
        # Assign the amplitude from the original state vector
        psi_new[j+1] = psi[i+1]
    end
    return psi_new
end