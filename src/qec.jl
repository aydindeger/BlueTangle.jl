"""
get_codeword(codestate::AbstractVectorS)
"""
get_codeword(codestate::AbstractVectorS) = fock_basis.(findall(x -> !isapprox(x, 0; atol=1000eps()), codestate) .- 1, get_N(codestate))

function stabilizers_to_generator(stabilizers::Vector)
    isempty(stabilizers) && throw(ArgumentError("at least one stabilizer is required"))
    rows = [strip.(split(strip(s), ",")) for s in stabilizers]
    num_qubits = length(first(rows))
    num_stabilizers = length(stabilizers)

    generator_matrix = zeros(Int, num_stabilizers, 2 * num_qubits)

    for (i, row) in enumerate(rows)
        length(row) == num_qubits || throw(ArgumentError("stabilizers must have the same number of qubits"))
        for (j, op) in enumerate(row)
            if op == "X"
                generator_matrix[i, j] = 1
            elseif op == "Z"
                generator_matrix[i, j+num_qubits] = 1
            elseif op == "Y"
                generator_matrix[i, j] = 1  # X part
                generator_matrix[i, j+num_qubits] = 1  # Z part
            elseif op != "I"
                throw(ArgumentError("stabilizers must contain only I, X, Y, or Z; got `$op`"))
            end
        end
    end

    return generator_matrix
end


"""
get_standard_form(G::AbstractMatrix)
"""
function get_standard_form(G::AbstractMatrix)
    iseven(size(G, 2)) || throw(ArgumentError("a binary stabilizer matrix must have 2n columns"))
    all(x -> x == 0 || x == 1, G) || throw(ArgumentError("the stabilizer matrix must be binary"))
    G = Int.(G)
    n, m = size(G, 2) ÷ 2, size(G, 1)
    m <= n || throw(ArgumentError("there cannot be more independent stabilizers than qubits"))
    all(iszero, mod.(G[:, 1:n] * G[:, n+1:end]' + G[:, n+1:end] * G[:, 1:n]', 2)) ||
        throw(ArgumentError("stabilizers must commute"))
    reduced_row_echelon(copy(G))[2] == m || throw(ArgumentError("stabilizers must be independent"))
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
        G_row_2 = hcat(zeros(Int, n - k - r, n), hcat(D, E_rref))
        G = vcat(G_row_1, G_row_2)
    end

    permute_matrix_E = Matrix{Int64}(la.I, n, n)
    len_E = size(E_transform_cols, 1)
    permute_matrix_E[n-len_E+1:n, n-len_E+1:n] .= E_transform_cols
    permute_matrix_inv = permute_matrix_E' * generator_transform_cols'

    if r != rank_of_rref(G)
        throw("something's wrong with the rank or standard form")
    end

    permute_vector_inv = [findfirst(permute_matrix_inv'[:, i] .== 1) for i in 1:n] #' is correct
    logical_indices = [permute_vector_inv[n-k+i] for i = 1:k]

    return G, permute_matrix_inv, permute_vector_inv, logical_indices
end

"""
get_XZ_logicals(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix)
"""
function get_XZ_logicals(G_standard::AbstractMatrix, permute_matrix_inv::AbstractMatrix)
    get_XZ_logicals!(G_standard::AbstractMatrix, permute_matrix_inv::AbstractMatrix, Dict())
end

"""
get_XZ_logicals!(G_standard::AbstractMatrix,permute_matrix_inv::AbstractMatrix,logical_XZ_ops::Dict)
"""
function get_XZ_logicals!(G_standard::AbstractMatrix, permute_matrix_inv::AbstractMatrix, logical_XZ_ops::Dict)

    r = rank_of_rref(G_standard)
    n, m = size(G_standard, 2) ÷ 2, size(G_standard, 1)
    k = n - m

    id, A1, A2, B, C1, C2, D, E = _extract_blocks(G_standard, r)

    logical_XZ_vecs = Dict()

    U0 = zeros(k, r)
    U2 = E'
    U3 = la.I(k)
    V1 = (E' * C1') + C2'
    U1 = zeros(k, m - r)
    V3 = zeros(k, k)
    #size: r m-r k
    X_vecs = [U0 U2 U3 V1 U1 V3]
    Z_vecs = [zeros(k, n) A2' U1 la.I(k)]

    # #permute
    ip = permute_matrix_inv#reduced_row_echelon_inverse(permute_matrix)
    Xvecs = Int.(mod.(hcat(X_vecs[:, 1:n] * ip, X_vecs[:, n+1:end] * ip), 2))
    Zvecs = Int.(mod.(hcat(Z_vecs[:, 1:n] * ip, Z_vecs[:, n+1:end] * ip), 2))

    for ki = 1:k
        logical_XZ_vecs["Z", ki] = sa.sparse(Zvecs[ki, :])
        logical_XZ_vecs["X", ki] = sa.sparse(Xvecs[ki, :])
        logical_XZ_ops["Z", ki] = logical_vec_2_ops(Zvecs[ki, :])
        logical_XZ_ops["X", ki] = logical_vec_2_ops(Xvecs[ki, :])
    end

    return logical_XZ_ops, sa.sparse(Xvecs), sa.sparse(Zvecs)

end

function logical_vec_2_ops(logical_vec::Vector)

    n = length(logical_vec) ÷ 2

    ops = []
    for i in 1:n
        x, z = mod(logical_vec[i], 2), mod(logical_vec[i+n], 2)
        if x == 1 || z == 1
            push!(ops, (x == 1 ? (z == 1 ? "Y" : "X") : "Z", i))
        end
    end

    return ops
end

function _extract_blocks(G_standard_form, r::Int) #r=rank
    #block sizes: r, m-r, k
    n, m = size(G_standard_form, 2) ÷ 2, size(G_standard_form, 1)
    k = n - m

    id = G_standard_form[1:r, 1:r]
    A1 = G_standard_form[1:r, r+1:n-k]
    A2 = G_standard_form[1:r, n-k+1:n]
    B = G_standard_form[1:r, n+1:n+r]
    C1 = G_standard_form[1:r, n+r+1:2n-k]
    C2 = G_standard_form[1:r, 2n-k+1:2n]

    D = G_standard_form[r+1:m, n+1:n+r]
    id2 = G_standard_form[r+1:m, n+r+1:2n-k]
    E = G_standard_form[r+1:m, 2n-k+1:2n]

    if (la.one(id) ≈ id) && (la.one(id2) ≈ id2)
        return id, A1, A2, B, C1, C2, D, E
    else
        throw("something's wrong with the standard form")
    end

end

"""
encoding_circuit_from_generator(generator_standard::AbstractMatrix,Xvecs::AbstractMatrixS,permute_vector_inv::AbstractVector)
"""
function encoding_circuit_from_generator(generator_standard::AbstractMatrix, Xvecs::AbstractMatrixS, permute_vector_inv::AbstractVector)

    m, n = size(generator_standard, 1), size(generator_standard, 2) ÷ 2
    r = rank_of_rref(generator_standard)
    k = n - m

    encoding_circuit = Vector{Op}()

    # First loop: Add CX gates based on logical_xs matrix
    for i in 1:k #goes over whole

        control_qubit = permute_vector_inv[n-k+i]

        # if Xvecs[i,control_qubit]==1
        #    #control qubit 
        # end

        for j in 1:n
            if control_qubit != j
                if Xvecs[i, j] == 1
                    push!(encoding_circuit, Op("CX", control_qubit, j))
                end
            end
        end

    end

    for i in 1:r

        control_qubit = permute_vector_inv[i]
        if generator_standard[i, i+n] == 0
            push!(encoding_circuit, Op("H", control_qubit))
        else
            push!(encoding_circuit, Op("H", control_qubit))
            push!(encoding_circuit, Op("S", control_qubit))
        end

        for j in 1:n
            if i != j
                if (generator_standard[i, j] == 1) && (generator_standard[i, j+n] == 0)
                    push!(encoding_circuit, Op("CX", control_qubit, permute_vector_inv[j]))
                elseif (generator_standard[i, j] == 0) && (generator_standard[i, j+n] == 1)
                    push!(encoding_circuit, Op("CZ", control_qubit, permute_vector_inv[j]))
                elseif (generator_standard[i, j] == 1) && (generator_standard[i, j+n] == 1)
                    push!(encoding_circuit, Op("CY", control_qubit, permute_vector_inv[j]))
                end
            end
        end
    end

    return encoding_circuit
end

"""
    get_ops_syndrome(generator::AbstractMatrix)

Measure each +1 Pauli generator with one ancilla appended after the data qubits.
The measurement bit is zero for eigenvalue +1 and one for eigenvalue -1.
"""
function get_ops_syndrome(generator_standard::AbstractMatrix)

    m, n = size(generator_standard, 1), size(generator_standard, 2) ÷ 2

    standard_generators_x = generator_standard[:, 1:n]
    standard_generators_z = generator_standard[:, n+1:end]

    syndrome_circuit = Vector{Op}()

    for i in 1:m
        push!(syndrome_circuit, Op("H", n + i))
    end

    for i in 1:m
        for j in 1:n
            if (standard_generators_x[i, j] == 1) && (standard_generators_z[i, j] == 1)
                push!(syndrome_circuit, Op("CY", n + i, j))
            elseif standard_generators_x[i, j] == 1
                push!(syndrome_circuit, Op("CX", n + i, j))
            elseif standard_generators_z[i, j] == 1
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

    Return a single-qubit Pauli with this syndrome, or an identity for zero
    syndrome. This decoder assumes at most one correctable Pauli error; a
    multi-qubit error can share the same syndrome and is not certified corrected.
"""
function identify_error(syndrome::AbstractVector, generator_standard::AbstractMatrix)

    if length(syndrome) != size(generator_standard, 1)
        throw(DimensionMismatch("syndrome length: expected $(size(generator_standard, 1)), got $(length(syndrome))"))
    end
    all(x -> x == 0 || x == 1, syndrome) || throw(ArgumentError("syndrome values must be zero or one"))
    all(iszero, syndrome) && return Op("I", 1; noisy=false)

    n = size(generator_standard, 2) ÷ 2
    # Iterate through all possible single qubit errors
    for qubit in 1:n
        for error_op in ("X", "Z", "Y")
            x = error_op in ("X", "Y")
            z = error_op in ("Z", "Y")
            # A single Pauli touches just two check-matrix columns.
            if all(i -> mod(x * generator_standard[i, n+qubit] + z * generator_standard[i, qubit], 2) == syndrome[i], eachindex(syndrome))
                return Op(error_op, qubit)
            end
        end
    end
    return "No single-qubit error matches this syndrome"
end

"""
physical_ops(name::String, logical_qubit::Int, logicals::Dict)
"""
physical_ops(name::String, logical_qubit::Int, logicals::Dict) = _physical_op_list(logicals[name, logical_qubit])

"""
physical_ops(name::String, logical_qubit::Int, target_qubit::Int, logicals::Dict)
"""
physical_ops(name::String, logical_qubit::Int, target_qubit::Int, logicals::Dict) = _physical_op_list(logicals[name, logical_qubit, target_qubit])

_physical_op_list(op::QuantumOps) = QuantumOps[op]
_physical_op_list(ops::AbstractVector) = QuantumOps[x isa QuantumOps ? x : Op(x...) for x in ops]

"""
physical_ops(op::QuantumOps, logicals::Dict)
"""
function physical_ops(op::QuantumOps, logicals::Dict)
    if isa(op, OpF)
        vals = logicals[op]
    elseif op.q == 1
        vals = logicals[op.name, op.qubit]
    elseif op.q == 2
        vals = logicals[op.name, op.qubit, op.target_qubit]
    end

    return vals isa QuantumOps ? vals : _physical_op_list(vals)
end

"""
physical_ops(ops::Vector{QuantumOps}, logicals::Dict)
"""
function physical_ops(ops::Vector{<:QuantumOps}, logicals::Dict)
    list_of_ops = Vector{QuantumOps}()
    for op = ops
        mapped = physical_ops(op, logicals)
        mapped isa QuantumOps ? push!(list_of_ops, mapped) : append!(list_of_ops, mapped)
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

# Conjugate a Hermitian Pauli P(x,z) by the Clifford gates used by the encoder.
# The returned sign is essential: binary row reduction alone forgets -P phases.
function _qec_conjugate_gate!(v, name, q, t=-1)
    n = length(v) ÷ 2
    x, z = v[q], v[n+q]
    if name == "H"
        v[q], v[n+q] = z, x
        return x * z == 1 ? -1 : 1
    elseif name == "S" || name == "SD"
        v[n+q] = xor(z, x)
        return x * (name == "S" ? z : 1-z) == 1 ? -1 : 1
    elseif name == "CX" || name == "CNOT"
        xt, zt = v[t], v[n+t]
        sign = x * zt * xor(xor(xt, z), 1) == 1 ? -1 : 1
        v[t], v[n+q] = xor(xt, x), xor(z, zt)
        return sign
    elseif name == "CZ"
        return _qec_conjugate_gate!(v, "H", t) *
               _qec_conjugate_gate!(v, "CX", q, t) *
               _qec_conjugate_gate!(v, "H", t)
    elseif name == "CY"
        return _qec_conjugate_gate!(v, "SD", t) *
               _qec_conjugate_gate!(v, "CX", q, t) *
               _qec_conjugate_gate!(v, "S", t)
    elseif name in ("I", "X", "Y", "Z")
        anticommutes = name == "X" ? z : name == "Z" ? x : name == "Y" ? xor(x, z) : 0
        return anticommutes == 1 ? -1 : 1
    end
    throw(ArgumentError("unsupported encoding Clifford gate: $name"))
end

function _qec_conjugate_pauli(v, ops; inverse=false)
    result = Int.(v)
    sign = 1
    for op in (inverse ? Iterators.reverse(ops) : ops)
        name = inverse && op.name == "S" ? "SD" : op.name
        sign *= _qec_conjugate_gate!(result, name, op.qubit, op.target_qubit)
    end
    return result, sign
end

function _qec_fix_encoding_signs(ops, generator, logical_indices)
    m, n = size(generator, 1), size(generator, 2) ÷ 2
    ancillas = setdiff(1:n, logical_indices)
    checks = zeros(Int, m, m)
    signs = zeros(Int, m)
    for row in 1:m
        decoded, sign = _qec_conjugate_pauli(generator[row, :], ops; inverse=true)
        all(iszero, decoded[1:n]) && all(iszero, decoded[n .+ logical_indices]) ||
            error("encoding circuit does not map the stabilizers to ancilla Z operators")
        checks[row, :] = decoded[n .+ ancillas]
        signs[row] = sign == -1
    end
    # Solve checks * bits = signs over GF(2), choosing the required input sector.
    _, rank, rows, cols = reduced_row_echelon(copy(checks))
    rank == m || error("encoding circuit has dependent ancilla checks")
    bits = mod.(cols * rows * signs, 2)
    preparation = [Op("X", ancillas[i]) for i in 1:m if bits[i] == 1]
    return vcat(preparation, ops)
end

_qec_append_ancillas(state::AbstractVector, m) = kron(state, zero_state(m))
function _qec_append_ancillas(rho::AbstractMatrix, m)
    ancillas = zero_state(m)
    return kron(rho, ancillas * ancillas')
end
function _qec_append_ancillas(psi::it.MPS, m)
    return it.MPS(vcat(collect(deepcopy(psi)), collect(zero_state(N_MPS(m)))))
end

"""
    StabilizerCode(stabilizers::Vector,logicals::Dict)

    Construct a code from independent commuting +1 Pauli generators. Rows may
    contain I, X, Y, and Z. Automatically derived logical gates respect the
    encoder's phases; explicitly supplied mappings are retained.

    `code.syndrome(state)` returns the measured state with reset ancillas, the
    syndrome bits in the input-generator order, and a candidate single-qubit
    correction. `code.correct` uses this bounded decoder. Pure-state `code.decode` requires
    the discarded qubits to be in |0⟩ after inverse encoding; it does not silently
    project an uncorrected state into a pure logical vector.

    Density matrices encode and decode by linear evolution and partial trace.
    `code.syndrome(rho)` samples a selective record, while `code.correct(rho)`
    sums the recovered syndrome sectors deterministically. MPS preparation,
    logical gates, syndrome extraction and correction stay in tensor form.

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
    ops_encoding::Vector #encoding circuit
    ops_syndrome::Vector #syndrome circuit
    info::NamedTuple
    ops::Function
    apply::Function
    encode::Function
    decode::Function
    syndrome::Function
    correct::Function

    function StabilizerCode(stabilizers::Union{Vector{String},String}, logicals::Dict=Dict(); d::Int=-1)

        if isa(stabilizers, String)
            stabilizers = [join(collect(line), ",") for line in split(chomp(stabilizers), '\n')]
        end

        generator = stabilizers_to_generator(stabilizers)
        n = size(generator, 2) ÷ 2
        m = length(stabilizers)
        k = n - m

        generator_standard, permute_matrix_inv, permute_vector_inv, logical_indices = get_standard_form(generator)
        _, Xvecs, Zvecs = get_XZ_logicals(generator_standard, permute_matrix_inv)
        ops_encoding = encoding_circuit_from_generator(generator_standard, Xvecs, permute_vector_inv)
        ops_encoding = _qec_fix_encoding_signs(ops_encoding, generator, logical_indices)
        ops_syndrome = get_ops_syndrome(generator)

        # Derive logical representatives from this encoder, retaining their signs.
        # Explicit user-supplied logical gates take precedence.
        logicals = Dict{Any,Any}(logicals)
        for i in 1:k, name in ("I", "X", "Y", "Z")
            v = zeros(Int, 2n)
            v[logical_indices[i]] = name in ("X", "Y")
            v[n+logical_indices[i]] = name in ("Z", "Y")
            encoded, sign = _qec_conjugate_pauli(v, ops_encoding)
            rep = logical_vec_2_ops(encoded)
            sign == -1 && pushfirst!(rep, ("phase", -gate.I, 1))
            get!(logicals, (name, i), rep)
            name == "X" && (Xvecs[i, :] = encoded)
            name == "Z" && (Zvecs[i, :] = encoded)
        end

        new_ops(name::String, logical_qubit::Int) = physical_ops(name, logical_qubit, logicals)
        new_ops(name::String, logical_qubit::Int, target_qubit::Int) = physical_ops(name, logical_qubit, target_qubit, logicals)
        new_ops(op::QuantumOps) = physical_ops(op, logicals)
        new_ops(ops::Vector{<:QuantumOps}) = physical_ops(ops, logicals)

        new_encode(state_init_sym::Symbol=:zero; random_op_count::Int=20, noise::Union{Bool,NoiseModel}=false, encoding::Vector=[]) = new_encode(n, state_init_sym; random_op_count=random_op_count, noise=noise, encoding=encoding)

        #mps
        function new_encode(N::Union{Int,Vector{<:its.Index}}, state_init_sym::Symbol=:zero; random_op_count::Int=20, noise::Union{Bool,NoiseModel}=false, encoding::Vector=[])

            (N isa Int ? N : length(N)) == n || throw(DimensionMismatch("encoding requires $n physical qubits"))

            state = qec_state_prep(N, logical_indices, state_init_sym; random_op_count=random_op_count)[1]

            ops = isempty(encoding) ? ops_encoding : encoding
            for o = ops
                state = apply(state, o; noise=noise)
            end

            return state

        end

        function new_encode(logical_state::AbstractVectorS; noise::Union{Bool,NoiseModel}=false, encoding::Vector=[])
            get_N(logical_state) == k || throw(DimensionMismatch("encoding requires a $k-qubit logical state"))
            state = sa.spzeros(promote_type(ComplexF64, eltype(logical_state)), 2^n)
            for basis in eachindex(logical_state)
                iszero(logical_state[basis]) && continue
                physical_basis = 0
                for i in 1:k
                    physical_basis |= (((basis-1) >> (k-i)) & 1) << (n-logical_indices[i])
                end
                state[physical_basis+1] = logical_state[basis]
            end
            return apply(isempty(encoding) ? ops_encoding : encoding, state; noise=noise)
        end

        function new_encode(logical_rho::AbstractMatrixS; noise::Union{Bool,NoiseModel}=false, encoding::Vector=[])
            size(logical_rho) == (2^k, 2^k) || throw(DimensionMismatch("encoding requires a $k-qubit logical density matrix"))
            indices = _partial_trace_offsets(fill(2, n), logical_indices) .+ 1
            rho = sa.spzeros(promote_type(ComplexF64, eltype(logical_rho)), 2^n, 2^n)
            rho[indices, indices] = logical_rho
            return apply(isempty(encoding) ? ops_encoding : encoding, rho; noise=noise)
        end

        function new_encode(logical_state::it.MPS; noise::Union{Bool,NoiseModel}=false, encoding::Vector=[])
            length(logical_state) == k || throw(DimensionMismatch("encoding requires a $k-qubit logical MPS"))
            state = _qec_append_ancillas(logical_state, n-k)
            state = it.movesites(state, collect(1:k) .=> logical_indices; cutoff=0.0)
            return apply(isempty(encoding) ? ops_encoding : encoding, state; noise=noise)
        end

        function new_decode(state::AbstractVectorS; noise::Union{Bool,NoiseModel}=false, decoding::Vector=[])

            get_N(state) in (n, n+m) || throw(DimensionMismatch("decoding requires $n data qubits, optionally followed by $m reset ancillas"))

            ops = isempty(decoding) ? ops_encoding' : decoding

            for o = ops
                state = apply(state, o; noise=noise)
            end

            return reduce_and_relabel_qubits(state, logical_indices)

        end

        function new_decode(rho::AbstractMatrixS; noise::Union{Bool,NoiseModel}=false, decoding::Vector=[])
            N = get_N(rho)
            N in (n, n+m) || throw(DimensionMismatch("decoding requires $n data qubits, optionally followed by $m ancillas"))
            rho = apply(isempty(decoding) ? ops_encoding' : decoding, rho; noise=noise)
            reduced = partial_trace(rho, fill(2, N), setdiff(1:N, logical_indices))
            permutation = [findfirst(==(q), sort(logical_indices)) for q in logical_indices]
            indices = _partial_trace_offsets(fill(2, k), permutation) .+ 1
            return reduced[indices, indices]
        end

        function new_decode(psi::it.MPS; noise::Union{Bool,NoiseModel}=false, decoding::Vector=[])
            length(psi) in (n, n+m) || throw(DimensionMismatch("decoding requires $n data qubits, optionally followed by $m ancillas"))
            psi = apply(isempty(decoding) ? ops_encoding' : decoding, psi; noise=noise)
            return reduce_and_relabel_qubits(psi, logical_indices)
        end

        function new_syndrome(state_encoded::Union{AbstractVectorS,AbstractMatrixS,it.MPS}; noise::Union{Bool,NoiseModel}=false)
            if get_N(state_encoded) != n
                throw(DimensionMismatch("syndrome extraction requires $n data qubits"))
            end

            state_encoded_ancillas = _qec_append_ancillas(state_encoded, m)
            state_encoded_ancillas, syndrome_vec = apply(ops_syndrome, state_encoded_ancillas;
                                                         noise=noise, track_measurements=true)
            op_error = identify_error(syndrome_vec, generator)
            for i in 1:m
                if syndrome_vec[i] == 1
                    state_encoded_ancillas = apply(state_encoded_ancillas, Op("X", n+i; noisy=false))
                end
            end

            return state_encoded_ancillas, syndrome_vec, op_error
        end

        function new_correct(state_encoded_ancillas::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, syndrome_vec::Vector; noise::Union{Bool,NoiseModel}=false)
            if 0 < d <= 2
                throw(ArgumentError("a code of distance $d cannot correct an arbitrary single-qubit error"))
            end
            get_N(state_encoded_ancillas) in (n, n+m) || throw(DimensionMismatch("correction requires $n data qubits, optionally followed by $m reset ancillas"))
            op_error = identify_error(syndrome_vec, generator)
            op_error isa Op || throw(ArgumentError(op_error))
            return all(iszero, syndrome_vec) ? state_encoded_ancillas : apply(state_encoded_ancillas, op_error; noise=noise)

        end

        function new_correct(state_encoded::Union{AbstractVectorS,it.MPS}; noise::Union{Bool,NoiseModel}=false)

            state_encoded_ancillas, syndrome_vec, _ = new_syndrome(state_encoded; noise=noise)
            return new_correct(state_encoded_ancillas, syndrome_vec; noise=noise)

        end

        function new_correct(rho::AbstractMatrixS; noise::Union{Bool,NoiseModel}=false)
            get_N(rho) == n || throw(DimensionMismatch("correction requires $n data qubits"))
            0 < d <= 2 && throw(ArgumentError("a code of distance $d cannot correct an arbitrary single-qubit error"))
            recovered = zero(sa.sparse(rho))
            if noise isa NoiseModel
                # Keep noisy ancilla gates coherent until their final measurement.
                joint = _qec_append_ancillas(sa.sparse(rho), m)
                joint = apply(ops_syndrome[1:end-m], joint; noise=noise)
                for outcome in 0:2^m-1
                    indices = (outcome+1):2^m:2^(n+m)
                    branch = joint[indices, indices]
                    real(la.tr(branch)) > 100eps() * abs(la.tr(rho)) || continue
                    recovered += new_correct(branch, int2bin(outcome, m); noise=noise)
                end
            else
                # Ideal syndrome extraction needs only data-space projectors.
                sectors = Tuple{AbstractMatrix,Vector{Int}}[(sa.sparse(rho), Int[])]
                identity = sa.sparse(1.0 * la.I, 2^n, 2^n)
                for stabilizer in stabilizers
                    S = string_to_matrix(stabilizer)
                    next_sectors = typeof(sectors)()
                    for (branch, bits) in sectors, bit in 0:1
                        projector = (identity + (1-2bit) * S) / 2
                        projected = projector * branch * projector'
                        real(la.tr(projected)) > 100eps() * abs(la.tr(rho)) || continue
                        push!(next_sectors, (projected, [bits; bit]))
                    end
                    sectors = next_sectors
                end
                for (branch, bits) in sectors
                    recovered += new_correct(branch, bits)
                end
            end
            return _qec_append_ancillas(recovered, m)
        end

        #one qubit
        function new_apply(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, name::String, qubit::Int; noise::Union{Bool,NoiseModel}=false)
            1 <= qubit <= k || throw(ArgumentError("logical qubit must be between 1 and $k"))
            ops = physical_ops(name, qubit, logicals)
            N = get_N(state)
            N in (n, n+m) || throw(DimensionMismatch("logical gates require an encoded $n-qubit state"))

            for o = ops
                state = apply(state, o; noise=noise)
            end

            return state
        end

        #two qubit
        function new_apply(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, name::String, logical_qubit::Int, target_qubit::Int; noise::Union{Bool,NoiseModel}=false)
            N = get_N(state)
            max_qt = max(logical_qubit, target_qubit)
            min_qt = min(logical_qubit, target_qubit)
            distance = abs(logical_qubit - target_qubit)

            if !(1 <= min_qt <= max_qt <= k)
                throw(ArgumentError("logical qubits must be between 1 and $k"))
            elseif distance == 0
                throw(ArgumentError("logical qubits must differ"))
            end
            N in (n, n+m) || throw(DimensionMismatch("logical gates require an encoded $n-qubit state"))

            keyexist = (name, logical_qubit, target_qubit) ∈ keys(logicals)

            if keyexist
                ops = physical_ops(name, logical_qubit, target_qubit, logicals)
            else
                throw("logical does not exist")
            end

            for o = ops
                state = apply(state, o; noise=noise)
            end

            return state
        end

        function new_apply(state::Union{AbstractVectorS,AbstractMatrixS,it.MPS}, op::Op; noise::Union{Bool,NoiseModel}=false)
            op.control == -2 || throw(ArgumentError("provide an explicit logical gate mapping for controlled operations"))
            return op.q == 1 ? new_apply(state, op.name, op.qubit; noise=noise) : new_apply(state, op.name, op.qubit, op.target_qubit; noise=noise)
        end

        info = (
            permute_vector=[findfirst(permute_matrix_inv[:, i] .== 1) for i in 1:n],
            len_logical=length(keys(logicals)),
            logical_XZ_vecs=(Xvecs, Zvecs)
        )

        return new(n, k, d, m, stabilizers, generator, generator_standard, permute_matrix_inv, permute_vector_inv, logical_indices, logicals, ops_encoding, ops_syndrome, info, new_ops, new_apply, new_encode, new_decode, new_syndrome, new_correct)

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


struct QECState
    logical::Union{AbstractVectorS,it.MPS}
    physical::Union{AbstractVectorS,it.MPS}
    encoded_state::Union{Nothing,AbstractVectorS,it.MPS}
    sym::Symbol
    is_mps::Bool
    dims::Vector{Int} #[dim_logical,dim_physical]
    ops_random::Vector{Op}
    qubit_mapping::Vector{Int}

    function QECState(
        n::Union{Int,Vector},
        logical_indices::Vector,
        state_init_sym::Union{Symbol,Vector}=:zero;
        random_op_count::Int=20,
        return_random::Bool=false,
        random_u3_count::Int=0,
        cutoff::Union{Nothing,Real}=nothing,
        encoding_circuit::Union{Nothing,AbstractVector{<:QuantumOps}}=nothing,
        qubit_mapping::AbstractVector{<:Integer}=Int[],
    )
        if return_random
            state_physical, state_logical, ops_random = qec_state_prep(n, logical_indices, state_init_sym; random_op_count=random_op_count, return_random=return_random, random_u3_count=random_u3_count)
        else
            state_physical, state_logical = qec_state_prep(n, logical_indices, state_init_sym; random_op_count=random_op_count, return_random=return_random, random_u3_count=random_u3_count)
            ops_random = Op[]
        end
        sym = isa(state_init_sym, Symbol) ? state_init_sym : :unknown
        is_mps = !isa(n, Int)

        if is_mps && cutoff !== nothing
            it.truncate!(state_logical; cutoff=cutoff)
            it.truncate!(state_physical; cutoff=cutoff)
        end

        final_qubit_mapping = Int.(collect(qubit_mapping))
        encoded_state = nothing
        if encoding_circuit !== nothing
            N = is_mps ? length(n) : n
            append!(final_qubit_mapping, length(final_qubit_mapping)+1:N)
            encoding_ops_no_swap, final_qubit_mapping = relabel_swap(collect(encoding_circuit), final_qubit_mapping)
            encoded_state_input = deepcopy(state_physical)
            if cutoff === nothing || !is_mps
                encoded_state = apply(encoding_ops_no_swap, encoded_state_input)
            else
                encoded_state = apply(encoding_ops_no_swap, encoded_state_input; cutoff=cutoff)
            end
        end

        dims = is_mps ? [it.maxlinkdim(state_logical), it.maxlinkdim(state_physical)] : [length(state_logical), length(state_physical)]

        new(state_logical, state_physical, encoded_state, sym, is_mps, dims, ops_random, final_qubit_mapping)
    end

end

"""
matchgate_ops(n::Int, depth::Int; include_z::Bool=true)

Create a nearest-neighbor matchgate (Gaussian) circuit. By default this is
number-conserving (GIVENS + RZ).
"""
function matchgate_ops(n::Int, depth::Int; include_z::Bool=true)
    ops = Op[]
    for layer in 1:depth
        if include_z
            for i in 1:n
                θ = randn() * π
                push!(ops, Op("RZ($(θ))", i))
            end
        end
        start = isodd(layer) ? 1 : 2
        for i in start:2:(n - 1)
            θ = randn() * π
            push!(ops, Op("GIVENS($(θ))", i, i + 1))
        end
    end
    return ops
end

"""
qec_state_prep(n::Union{Int,Vector},logical_indices::Vector,state_init_sym::Symbol=:zero;random_op_count::Int=20)

If you want `QECState` to also build the encoded state, pass `encoding_circuit=...`
to the constructor so it can relabel away SWAPs and return the final qubit mapping.
state_init_sym supports :zero, :one, :plus, :minus, :zero_plus, :random, :gaussian (or :matchgate)
:gaussian starts from a random half-filled computational basis state on the logical qubits.
random_u3_count randomly picks count logical qubits and prepares them as random U3|0⟩ product states.
"""
function qec_state_prep(
    n::Union{Int,Vector},
    logical_indices::Vector,
    state_init_sym::Union{Symbol,Vector}=:zero;
    random_op_count::Int=20,
    return_random::Bool=false,
    random_u3_count::Int=0,
)
    # state_init_sym=:random

    mps_bool = isa(n, Integer) ? false : true

    len_k = length(logical_indices)
    N = mps_bool ? length(n) : n
    all(i -> i isa Integer && 1 <= i <= N, logical_indices) ||
        throw(ArgumentError("logical indices must lie between 1 and $N"))
    length(unique(logical_indices)) == len_k || throw(ArgumentError("logical indices must be unique"))
    random_op_count >= 0 || throw(ArgumentError("random_op_count must be nonnegative"))
    if random_u3_count < 0 || random_u3_count > len_k
        throw(ArgumentError("random_u3_count must be between 0 and $(len_k)"))
    end
    if random_u3_count > 0 && state_init_sym in (:random, :gaussian, :matchgate)
        throw(ArgumentError("random_u3_count is only supported for product initial states"))
    end
    random_u3_sites = random_u3_count == 0 ? Int[] : sort(Int.(sb.sample(1:len_k, random_u3_count; replace=false)))

    k = mps_bool == false ? len_k : n[1:len_k]

    state = zero_state(n)
    state_logical = zero_state(k)
    ops_random_u3 = Op[]

    if state_init_sym isa Vector
        length(state_init_sym) == len_k || throw(DimensionMismatch("provide one bit per logical qubit"))
        all(x -> x == 0 || x == 1, state_init_sym) || throw(ArgumentError("a product-state vector must contain only zero and one"))
        for i in 1:len_k
            i in random_u3_sites && continue
            if state_init_sym[i] == 1
                state = Op("X", logical_indices[i]) * state
                state_logical = Op("X", i) * state_logical
            end
        end
    elseif (state_init_sym == :zeros) || (state_init_sym == :zero)
        nothing
    elseif (state_init_sym == :ones) || (state_init_sym == :one) #check this

        for i = 1:len_k
            i in random_u3_sites && continue
            state = Op("X", logical_indices[i]) * state
            state_logical = Op("X", i) * state_logical
        end

    elseif state_init_sym == :plus
        for i = 1:len_k
            i in random_u3_sites && continue
            state = Op("H", logical_indices[i]) * state
            state_logical = Op("H", i) * state_logical
        end
    elseif state_init_sym == :minus #check this

        for i = 1:len_k
            i in random_u3_sites && continue
            state = Op("X", logical_indices[i]) * state
            state = Op("H", logical_indices[i]) * state
            state_logical = Op("X", i) * state_logical
            state_logical = Op("H", i) * state_logical
        end

    elseif state_init_sym == :zero_plus
        # |0⟩⊗|+⟩^(k-1): first qubit in |0⟩, rest in |+⟩
        for i = 2:len_k
            i in random_u3_sites && continue
            state = Op("H", logical_indices[i]) * state
            state_logical = Op("H", i) * state_logical
        end

    elseif state_init_sym == :random

        ops_random = BlueTangle.random_ops(len_k, random_op_count)

        for o = ops_random
            state_logical = o * state_logical
            if o.q == 1
                state = Op(o.name, logical_indices[o.qubit]) * state
            elseif o.q == 2
                state = Op(o.name, logical_indices[o.qubit], logical_indices[o.target_qubit]) * state
            end
        end
        if return_random == true
            return state, state_logical, ops_random
        end
    elseif state_init_sym == :gaussian || state_init_sym == :matchgate
        n_filled = len_k ÷ 2
        if isodd(len_k) && rand(Bool)
            n_filled += 1
        end
        occ = n_filled == 0 ? Int[] : sb.sample(1:len_k, n_filled; replace=false)
        for i in occ
            state = Op("X", logical_indices[i]) * state
            state_logical = Op("X", i) * state_logical
        end

        ops_random = matchgate_ops(len_k, random_op_count)

        for o = ops_random
            state_logical = o * state_logical
            if o.q == 1
                state = Op(o.name, logical_indices[o.qubit]) * state
            elseif o.q == 2
                state = Op(o.name, logical_indices[o.qubit], logical_indices[o.target_qubit]) * state
            end
        end
        if return_random == true
            return state, state_logical, ops_random
        end
    else
        throw(ArgumentError("invalid state name"))
    end

    for i in random_u3_sites
        theta = acos(1 - 2 * rand())
        phi = 2 * pi * rand()
        op = Op("U3($(theta),$(phi),0)", i)
        push!(ops_random_u3, op)
        state_logical = op * state_logical
        state = Op(op.name, logical_indices[i]) * state
    end

    if return_random == true
        return state, state_logical, ops_random_u3
    end

    return state, state_logical
end

"""
reduce_and_relabel_qubits(psi::AbstractVectorS, logical_indices::Vector,qubit_mapping::Vector=Vector[])

Keep the listed qubits in the requested order. Every discarded qubit must be
in |0⟩; use `partial_trace` when discarding an occupied or entangled subsystem.
"""
function reduce_and_relabel_qubits(psi::AbstractVectorS, logical_indices::Vector, qubit_mapping::Vector=Vector[])
    n = get_N(psi)
    k = length(logical_indices)          # Number of qubits to keep
    psi_new = sa.spzeros(ComplexF64, 2^k)   # New state vector for k qubits

    if !isempty(qubit_mapping)
        qubit_mapping = [qubit_mapping; length(qubit_mapping)+1:n]
        sort(qubit_mapping) == collect(1:n) || throw(ArgumentError("qubit_mapping must permute its qubits without repetition"))
        logical_indices = [findfirst(x -> x == y, qubit_mapping) for y = logical_indices]
    end
    all(i -> i isa Integer && 1 <= i <= n, logical_indices) || throw(ArgumentError("invalid retained qubit index"))
    length(unique(logical_indices)) == k || throw(ArgumentError("retained qubit indices must be unique"))

    # Iterate over all basis states in the reduced Hilbert space
    for j in 0:(2^k-1)
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
    isapprox(sum(abs2, psi_new), sum(abs2, psi); atol=1e-12, rtol=1e-10) ||
        throw(ArgumentError("discarded qubits are not in |0⟩; correct/reset them first, or use partial_trace"))
    return psi_new
end

function reduce_and_relabel_qubits(psi::it.MPS, logical_indices::Vector, qubit_mapping::Vector=Int[])
    n, k = length(psi), length(logical_indices)
    if !isempty(qubit_mapping)
        mapping = [qubit_mapping; length(qubit_mapping)+1:n]
        sort(mapping) == collect(1:n) || throw(ArgumentError("qubit_mapping must permute its qubits without repetition"))
        logical_indices = [findfirst(==(q), mapping) for q in logical_indices]
    end
    all(q -> q isa Integer && 1 <= q <= n, logical_indices) && allunique(logical_indices) ||
        throw(ArgumentError("retained qubits must be distinct valid indices"))
    order = [logical_indices; setdiff(1:n, logical_indices)]
    reordered = it.movesites(psi, order .=> collect(1:n); cutoff=0.0)
    sites = get_M(reordered)
    environment = its.ITensor(1.0)
    for q in n:-1:k+1
        environment *= reordered[q] * its.state(sites[q], "0")
    end
    if k == 0
        result = ComplexF64[its.scalar(environment)]
    else
        tensors = [reordered[q] for q in 1:k]
        tensors[end] *= environment
        result = it.MPS(tensors)
    end
    isapprox(la.norm(result), la.norm(psi); atol=1e-12, rtol=1e-10) ||
        throw(ArgumentError("discarded qubits are not in |0⟩; correct/reset them first"))
    return result
end
