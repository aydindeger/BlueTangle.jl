"""
_get_codeword
"""
get_codeword(codestate::sa.SparseVector)=fock_basis.(findall(x->!isapprox(x,0;atol=1000eps()),codestate).-1,get_N(codestate))


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

    function StabilizerCode(n::Int,k::Int,d::Int,stabilizers::Vector,logicals::Dict)

        if length(first(logicals)) != k
            throw("k and logicals are incompatible")
        end

        ##find init from stabilizers
        stabilizer_matrix=string_to_matrix.(stabilizers)
        e,v=la.eigen(Matrix(mean(stabilizer_matrix)))
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
        len_stabilizers=length(stabilizers),
        len_codestates=len_codestates
        )

        return new(n,k,d,stabilizers,logicals,codestates,codewords,stabilizer_matrix,stats,new_ops,new_apply,new_encode,new_decode,new_syndrome,new_correct,new_circuit);

    end

end