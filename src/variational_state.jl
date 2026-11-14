# Apply only the local matrix, without constructing a 2^N × 2^N operator.
# Qubit 1 is the most significant state bit; the matrix operand order is
# (op.qubit, op.target_qubit), even when those sites are reversed or nonlocal.
function _variational_gate(state::AbstractVector, op::Op, pars=())
    matrix = _variational_gate_matrix(op, pars)
    result = Vector{promote_type(eltype(state), eltype(matrix))}(state)
    return _variational_gate!(result, matrix, op.qubit, op.target_qubit, op.control)
end

# The circuit evaluator owns this vector. Reuse it unless a parameterized gate
# introduces a wider scalar type (for example complex ForwardDiff amplitudes).
function _variational_gate!(state::Vector, op::Op, pars=())
    matrix = _variational_gate_matrix(op, pars)
    T = promote_type(eltype(state), eltype(matrix))
    result = T === eltype(state) ? state : Vector{T}(state)
    return _variational_gate!(result, matrix, op.qubit, op.target_qubit, op.control)
end

function _variational_gate_matrix(op::Op, pars)
    op.type == "🔬" && throw(ArgumentError("the local variational kernel does not sample measurements"))
    return op.mat isa Function ? op.mat(pars...) : op.mat
end

function _variational_gate!(state::Vector, matrix::AbstractMatrix,
                            qubit::Int, target::Int, control::Int)
    dimension = length(state)
    ispow2(dimension) || throw(DimensionMismatch("state length must be a positive power of two"))
    N = trailing_zeros(dimension)
    1 <= qubit <= N || throw(ArgumentError("qubit must be within the state register"))
    target == -1 || (1 <= target <= N && target != qubit) ||
        throw(ArgumentError("target must be a distinct qubit within the state register"))
    control == -2 || (1 <= control <= N && control != qubit && control != target) ||
        throw(ArgumentError("control must be a distinct qubit within the state register"))
    local_dimension = target == -1 ? 2 : 4
    size(matrix) == (local_dimension, local_dimension) ||
        throw(DimensionMismatch("gate matrix dimensions do not match its operands"))

    first_bit = 1 << (N - qubit)
    second_bit = target == -1 ? 0 : 1 << (N - target)
    control_bit = control == -2 ? 0 : 1 << (N - control)
    if second_bit == 0
        @inbounds for block in 0:2first_bit:length(state)-1, offset in 0:first_bit-1
            base = block + offset
            base & control_bit == control_bit || continue
            i, j = base + 1, base + first_bit + 1
            a, b = state[i], state[j]
            state[i] = matrix[1,1]*a + matrix[1,2]*b
            state[j] = matrix[2,1]*a + matrix[2,2]*b
        end
    else
        operated_bits = first_bit | second_bit
        @inbounds for base in 0:length(state)-1
            iszero(base & operated_bits) && base & control_bit == control_bit || continue
            i = base + 1
            j = base + second_bit + 1
            k = base + first_bit + 1
            l = base + operated_bits + 1
            a, b, c, d = state[i], state[j], state[k], state[l]
            state[i] = matrix[1,1]*a + matrix[1,2]*b + matrix[1,3]*c + matrix[1,4]*d
            state[j] = matrix[2,1]*a + matrix[2,2]*b + matrix[2,3]*c + matrix[2,4]*d
            state[k] = matrix[3,1]*a + matrix[3,2]*b + matrix[3,3]*c + matrix[3,4]*d
            state[l] = matrix[4,1]*a + matrix[4,2]*b + matrix[4,3]*c + matrix[4,4]*d
        end
    end
    return state
end
