# A small adapter gives BlueTangle methods a concrete dispatch target and keeps
# its most-significant-bit convention separate from PhasePoly's native ordering.
struct PhasePolyState{S}
    sim::S
    phase::ComplexF64
end

function phase_gates!(program, angle, qubits)
    step = length(qubits) == 1 ? pi/8 : pi/4
    units = round(Int, angle/step)
    isapprox(angle, units*step; atol=1e-10, rtol=0) ||
        throw(ArgumentError("PhasePoly requires single-qubit phases in multiples of π/8 and controlled phases in multiples of π/4."))
    names = length(qubits) == 1 ? ("Z", "S", "T", "rT") : ("CZ", "CS", "CT")
    units = mod(units, 1 << length(names))
    for (i, name) in enumerate(names)
        units & (1 << (length(names)-i)) != 0 && push!(program, (name, qubits))
    end
end

"""
    phasepoly(ops; nqubits=get_N(ops), kwargs...)

Convert a noiseless supported circuit to PhasePoly's compact representation.
`ops` may be a vector of `QuantumOps`, tuple shorthands such as
`[("H", 1), ("CX", 1, 2)]`, or a mix of both. Load `PhasePoly` first to enable
this optional extension. Physical qubit labels are unchanged; `to_state` and
`amplitude` retain BlueTangle's basis ordering. The returned adapter exposes the
native simulator as `.sim`.

Supported gates are H, affine permutations, Y, diagonal Clifford/T gates, and
P/U1/RZ phases in multiples of π/8 (controlled phases in multiples of π/4).
Noise, measurement, unbound parameters, and unsupported gates raise errors.
Leading and trailing H layers use PhasePoly's efficient paths; interior H gates
require 2^h branches, with PhasePoly's `max_h=12` default. Level-4 phases can also
require exponential branching. Other PhasePoly cost guards pass through `kwargs`.
Pauli expectations use the native evaluator; `to_state` explicitly materializes
a sparse vector subject to PhasePoly's `max_nonzeros` limit.
"""
function BlueTangle.phasepoly(ops::Vector{<:QuantumOps}; N::Int=get_N(ops),
                              nqubits::Union{Nothing,Int}=nothing, kwargs...)
    nqubits = isnothing(nqubits) ? N : nqubits
    nqubits > 0 && nqubits >= get_N(ops) || throw(ArgumentError("nqubits must cover the circuit."))
    program = Tuple{String,Vector{Int}}[]
    global_phase = 1.0 + 0.0im
    aliases = Dict("SD"=>"Sdg", "TD"=>"Tdg", "CSD"=>"CSdg", "CTD"=>"CTdg", "CNOT"=>"CX")
    for op in ops
        op isa Op && op.type != "🔬" && op.mat isa AbstractMatrix ||
            throw(ArgumentError("PhasePoly accepts only bound unitary gate operations."))
        op.mat ≈ gates(op.name) || throw(ArgumentError("Custom matrix for $(op.name) cannot be translated by name."))
        name = BlueTangle._clean_name(op.name)
        q, t, c = op.qubit, op.target_qubit, op.control
        if name in ("P", "U1", "RZ", "CP")
            angle = only(BlueTangle._get_name_and_arguments(op.name)[3])
            if name == "CP"
                c < 0 || throw(ArgumentError("A further-controlled CP gate is unsupported."))
                phase_gates!(program, angle, [q,t])
            elseif c > 0
                name == "RZ" && phase_gates!(program, -angle/2, [c])
                phase_gates!(program, angle, [c,q])
            else
                name == "RZ" && (global_phase *= cis(-angle/2))
                phase_gates!(program, angle, [q])
            end
            continue
        end
        if c > 0
            if op.q == 1 && name in ("I", "X", "Y", "Z", "S", "SD", "T", "TD")
                name == "I" && continue
                if name == "Y"
                    append!(program, [("Sdg",[q]), ("CX",[c,q]), ("S",[q])])
                else
                    controlled = Dict("X"=>"CX", "Z"=>"CZ", "S"=>"CS", "SD"=>"CSdg", "T"=>"CT", "TD"=>"CTdg")
                    push!(program, (controlled[name],[c,q]))
                end
            elseif op.q == 2 && name in ("CZ", "CS", "CSD")
                push!(program, (name == "CZ" ? "CCZ" : name == "CS" ? "CCS" : "CCSdg",[c,q,t]))
            else
                throw(ArgumentError("PhasePoly does not support controlled $(op.name)."))
            end
        elseif name in ("I", "CI")
            continue
        elseif name == "CY"
            append!(program, [("Sdg",[t]), ("CX",[q,t]), ("S",[t])])
        elseif name == "HSP"
            append!(program, [("Sdg",[q]), ("H",[q])])
        elseif name == "HY"
            append!(program, [("S",[q]), ("H",[q])])
        elseif name in ("H", "X", "Y", "Z", "S", "SD", "T", "TD", "CX", "CNOT", "SWAP", "CZ", "CS", "CSD", "CT", "CTD")
            push!(program, (get(aliases,name,name), op.q == 1 ? [q] : [q,t]))
        else
            throw(ArgumentError("PhasePoly does not support $(op.name)."))
        end
    end

    leading = 0
    while leading < length(program) && program[leading+1][1] == "H"
        leading += 1
    end
    trailing = length(program)
    final_h = Set{Int}()
    while trailing > leading && program[trailing][1] == "H"
        q = only(program[trailing][2])
        q in final_h ? delete!(final_h,q) : push!(final_h,q)
        trailing -= 1
    end
    native = PhasePoly.simulate(program[1:leading],program[leading+1:trailing];
                                n=nqubits, final_h=sort!(collect(final_h)), kwargs...)
    return PhasePolyState(native,global_phase)
end

BlueTangle.phasepoly(ops::Union{BlueTangle._OpList,AbstractVector{Any}}; kwargs...) =
    BlueTangle.phasepoly(BlueTangle._as_ops(ops); kwargs...)

BlueTangle.get_N(state::PhasePolyState) = state.sim.n

function BlueTangle.to_state(state::PhasePolyState; kwargs...)
    native = PhasePoly.statevector(state.sim; kwargs...)
    indices, values = SparseArrays.findnz(native)
    n = get_N(state)
    reordered = [Int(bitreverse(UInt(i-1)) >> (8sizeof(Int)-n))+1 for i in indices]
    return SparseArrays.sparsevec(reordered,state.phase .* values,length(native))
end

BlueTangle.expect(state::PhasePolyState, axis::String, qubit::Int) = PhasePoly.expect(state.sim,(axis,qubit))
function BlueTangle.expect(state::PhasePolyState, pauli::String)
    return length(pauli) == 1 ? [PhasePoly.expect(state.sim,(pauli,q)) for q in 1:get_N(state)] :
                               PhasePoly.expect(state.sim,pauli)
end
BlueTangle.expect(state::PhasePolyState, term::Tuple{Any,Integer}) = PhasePoly.expect(state.sim,term)
BlueTangle.expect(state::PhasePolyState, terms::AbstractVector{<:Tuple}) = PhasePoly.expect(state.sim,terms)
function BlueTangle.correlation(state::PhasePolyState, operators::String, qubits::Vector)
    names = strip.(split(operators,','))
    length(names) == length(qubits) || throw(DimensionMismatch("One Pauli operator is required per qubit."))
    return PhasePoly.expect(state.sim,collect(zip(names,qubits)))
end

function BlueTangle.amplitude(state::PhasePolyState, index::Integer; kwargs...)
    n = get_N(state)
    0 <= index < (big(1) << n) || throw(ArgumentError("Basis index is outside the state."))
    bits = [Int((index >> (n-q)) & 1) for q in 1:n]
    return state.phase * PhasePoly.amplitude(state.sim,bits; kwargs...)
end
BlueTangle.amplitude(state::PhasePolyState, bits::Union{AbstractVector,AbstractString}; kwargs...) =
    state.phase * PhasePoly.amplitude(state.sim,bits; kwargs...)

function BlueTangle.to_rho(state::PhasePolyState; kwargs...)
    vec = BlueTangle.to_state(state; kwargs...)
    return vec * vec'
end

BlueTangle.sample_exact(state::PhasePolyState; kwargs...) =
    BlueTangle.sample_exact(BlueTangle.to_state(state; kwargs...))
BlueTangle.sample(state::PhasePolyState, shots::Int=1; kwargs...) =
    BlueTangle.sample(BlueTangle.to_state(state; kwargs...), shots)
BlueTangle.sample_bit(state::PhasePolyState, shots::Int=1; basis="Z", kwargs...) =
    BlueTangle.sample_bit(BlueTangle.to_state(state; kwargs...), shots; basis=basis)
BlueTangle.sample_state(state::PhasePolyState, shots::Int; kwargs...) =
    BlueTangle.get_probs_from_sample(BlueTangle.sample(state, shots; kwargs...), BlueTangle.get_N(state))

function BlueTangle.measure(state::PhasePolyState, number_of_experiment::Int=-1;
                            label="state to measurement", kwargs...)
    number_of_experiment == -1 || number_of_experiment > 0 ||
        throw(ArgumentError("shots must be positive, or -1 for exact probabilities"))
    n = BlueTangle.get_N(state)
    bitstr, probabilities = number_of_experiment == -1 ?
        BlueTangle.sample_exact(state; kwargs...) :
        BlueTangle.sample_state(state, number_of_experiment; kwargs...)
    return BlueTangle._measurement_from_probs(bitstr, probabilities, n, number_of_experiment, "0", label)
end

BlueTangle.fidelity(a::PhasePolyState, b::AbstractVector; kwargs...) =
    BlueTangle.fidelity(BlueTangle.to_state(a; kwargs...), b)
BlueTangle.fidelity(a::AbstractVector, b::PhasePolyState; kwargs...) =
    BlueTangle.fidelity(a, BlueTangle.to_state(b; kwargs...))
BlueTangle.fidelity(a::PhasePolyState, b::PhasePolyState; kwargs...) =
    BlueTangle.fidelity(BlueTangle.to_state(a; kwargs...), BlueTangle.to_state(b; kwargs...))
