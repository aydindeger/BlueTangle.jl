# Density-matrix ensemble evolution.

"""
    apply(rho::AbstractMatrix, op::QuantumOps; noise=false)

Evolve a density matrix. By default, measurements and conditional operations
sum their unnormalized branches, and `MR` averages over X, Y and Z bases.
Use `run(rho, op)` to sample a measurement outcome and return the conditional
state and its record. Noise channels always apply the complete Kraus sum.
The input matrix is not modified.
"""
function apply(rho::AbstractMatrixS, op::QuantumOps;
               noise::Union{NoiseModel,Bool}=false, track_measurements::Bool=false, kwargs...)
    _unsupported_keyword_check(kwargs, [:pars])
    N = get_N(rho)
    bits = Int[]
    if op isa OpF
        if op.data isa Vector{<:QuantumOps}
            if track_measurements
                rho, bits = apply(rho, op.data; noise=noise, track_measurements=true, kwargs...)
            else
                rho = apply(rho, op.data; noise=noise, kwargs...)
            end
        else
            rho = op.apply(rho)
        end
    elseif op isa OpQC
        rho = op.apply(rho)
    elseif op.type == "🔬"
        if op isa ifOp
            rho, record = op.born_apply(rho, noise; track_measurements, kwargs...)
            track_measurements && append!(bits, record)
        else
            rho, bit = _measure_density(rho, op.name, op.qubit; selective=track_measurements)
            track_measurements && push!(bits, bit)
        end
    else
        operator = op.type == "f" ? op.expand(N, kwargs[:pars]...) : op.expand(N)
        rho = operator * rho * operator'
    end
    noise isa NoiseModel && (rho = apply_noise(rho, op, noise))
    return track_measurements ? (rho, bits) : rho
end


# Project in the physical measurement basis. Leaving branches unnormalized
# keeps both Born probabilities available for deterministic ensemble evolution.
function _density_measurement_branches(rho::AbstractMatrixS, name::String, qubit::Int)
    rotation = _measurement_mat(name)
    N = get_N(rho)
    return map((gate.P0, gate.P1)) do projector
        extended = hilbert(N, rotation' * projector * rotation, qubit)
        extended * rho * extended'
    end
end

function _measure_density(rho::AbstractMatrixS, name::String, qubit::Int; selective=false)
    if selective
        branches = _density_measurement_branches(rho, _resolve_measurement_name(name), qubit)
        probabilities = max.(0.0, real.(la.tr.(branches)))
        total = sum(probabilities)
        total > 0 || throw(ArgumentError("cannot sample a density matrix with zero trace"))
        bit = rand() < probabilities[1] / total ? 0 : 1
        return branches[bit + 1] / probabilities[bit + 1], bit
    end

    bases = uppercase(name) in ("MR", "M(R)") ? ("MX", "MY", "MZ") : (name,)
    result = zero(rho)
    for basis in bases
        result += sum(_density_measurement_branches(rho, basis, qubit)) / length(bases)
    end
    return result, nothing
end

_born_measure(rho::AbstractMatrixS, op::QuantumOps) = first(_measure_density(rho, op.name, op.qubit))

"""
    born_measure_Z(N, rho::AbstractMatrix, qubit)

Return the unconditioned Z-measurement channel, retaining both outcomes in a
density matrix. This overload returns only the matrix, with no sampled bit.
To sample an outcome and return a conditional matrix with its record, use
`run(rho, Op("MZ", qubit))`.
"""
born_measure_Z(N::Int, rho::AbstractMatrixS, qubit::Int) =
    first(_measure_density(rho, "MZ", qubit))


"""
    to_rho(state::AbstractVector)

Form the pure-state density matrix `state * state'`. To simulate the complete
noisy ensemble, apply operations to this density matrix:
`apply(to_rho(initial), ops; noise=noise)`.
"""
to_rho(state::AbstractVectorS) = state * state'

function sample_exact(rho::AbstractMatrixS)
    probabilities = real.(la.diag(rho))
    indices = findall(!iszero, probabilities)
    return indices .- 1, probabilities[indices]
end

function sample(rho::AbstractMatrixS, shots)
    outcomes, probabilities = sample_exact(rho)
    return sb.sample(outcomes, sb.Weights(probabilities), shots)
end
