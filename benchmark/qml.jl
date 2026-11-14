# Compare local statevector evolution with the previous expanded-gate approach.
# Run from the repository root: julia --project benchmark/qml.jl
# Both paths use the same circuits, initial states, gradients, and Adam updates.
# Compilation is warmed before timing; results depend on Julia and the machine.
using BlueTangle, ForwardDiff, LinearAlgebra, Statistics, Printf

sizes = (8, 10)
layers = 3
repetitions = 5
training_steps = 3

# This reference circuit contains only one-parameter RY/RZ gates and fixed CXs.
function expanded_apply(parameters, ops, initial)
    state = initial
    N = get_N(state)
    parameter = 1
    for op in ops
        if op.mat isa Function
            state = op.expand(N, parameters[parameter]) * state
            parameter += 1
        else
            state = op.expand(N) * state
        end
    end
    return state
end

function warmed_measurements(f, repetitions)
    f()
    seconds = median([(@elapsed f()) for _ in 1:repetitions])
    bytes = median([(@allocated f()) for _ in 1:repetitions])
    return seconds, bytes
end

function compare_qml(N, layers, repetitions, training_steps)
    ops = generate_ansatz_circuit(N, layers, ["RY", "RZ", "CX"])
    parameters = collect(range(0.17, step=0.013, length=2*N*layers))
    initial = Vector(zero_state(N))
    # Z expectation on q1; q1 is the most significant basis-state bit.
    z = [iszero((i-1) & (1 << (N-1))) ? 1.0 : -1.0 for i in eachindex(initial)]
    reference_state(p) = expanded_apply(p, ops, initial)
    local_state(p) = variational_apply(p, ops; init=initial)
    reference_loss(p) = sum(z .* abs2.(reference_state(p)))
    local_loss(p) = sum(z .* abs2.(local_state(p)))
    reference_gradient() = ForwardDiff.gradient(reference_loss, parameters)
    local_gradient() = ForwardDiff.gradient(local_loss, parameters)
    reference_training() = train(reference_loss, parameters; optimizer=:adam,
                                 iterations=training_steps, history=false)
    local_training() = train(local_loss, parameters; optimizer=:adam,
                             iterations=training_steps, history=false)

    @assert local_state(parameters) ≈ reference_state(parameters)
    @assert local_gradient() ≈ reference_gradient()
    old_result, new_result = reference_training(), local_training()
    @assert new_result.parameters ≈ old_result.parameters
    @assert new_result.loss ≈ old_result.loss

    println("\nN=$N, $layers layers, $(length(parameters)) parameters; median of $repetitions warmed calls")
    @printf("%-18s %12s %12s %10s %12s %12s\n", "Operation", "Expanded ms", "Local ms", "Speedup", "Expanded MiB", "Local MiB")
    cases = (("Forward", () -> reference_state(parameters), () -> local_state(parameters)),
             ("Gradient", reference_gradient, local_gradient),
             ("Adam × $training_steps", reference_training, local_training))
    for (label, reference, local_call) in cases
        old_time, old_bytes = warmed_measurements(reference, repetitions)
        new_time, new_bytes = warmed_measurements(local_call, repetitions)
        @printf("%-18s %12.3f %12.3f %9.1f× %12.3f %12.3f\n",
                label, 1000old_time, 1000new_time, old_time/new_time, old_bytes/2^20, new_bytes/2^20)
    end
end

for N in sizes
    compare_qml(N, layers, repetitions, training_steps)
end
