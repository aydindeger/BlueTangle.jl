# Regenerate the figures displayed in the documentation using the current API.
# Run from an environment with BlueTangle and PyPlot installed:
#   MPLBACKEND=Agg julia --project=/path/to/plot-env docs/generate_figures.jl
# Seeds make these small teaching examples reproducible; they are not research
# datasets. No hardware jobs are submitted. Outputs replace docs/src/assets/figs.
using BlueTangle, PyPlot, Random, LinearAlgebra, Statistics

figure_directory = joinpath(@__DIR__, "src", "assets", "figs")
function save_plot(name, input...; kwargs...)
    fig = plotq(input...; backend=:pyplot, show=false, kwargs...)
    savefigure(fig, joinpath(figure_directory, name); dpi=150)
    close(fig)
end

function save_scientific_plot(name, fig, ax)
    for side in ("top", "right")
        ax.spines[side].set_visible(false)
    end
    ax.grid(axis="y", color="#e5e9ef", linewidth=0.8)
    ax.set_axisbelow(true)
    fig.tight_layout(pad=1.1)
    savefigure(fig, joinpath(figure_directory, name); dpi=150)
    close(fig)
end

## Plotting gallery: feedback and exact noisy Bell-state probabilities.
feedback_ops = QuantumOps[
    Op("H", 1), Op("CX", 1, 2), Op("RY(0.37)", 3),
    Op("SWAP", 2, 3), Op("MZ", 2; if1="X"),
]
save_plot("circuit_feedback.png", feedback_ops)
bell_ops = [Op("H", 1), Op("CX", 1, 2)]
rho0 = zero_state(2) * zero_state(2)'
ideal = measure(apply(rho0, bell_ops); label="Ideal")
noisy = measure(apply(rho0, bell_ops; noise=NoiseModel("depolarizing", 0.08)); label="With gate noise")
save_plot("measurement_comparison.png", [ideal, noisy], ["Ideal", "With gate noise"])

## Examples 6 and 8: conditional measurement and a random Clifford circuit.
conditional_ops = QuantumOps[
    Op("X", 1), Op("MX", 1; if0="Z", if1="H"),
    Op("CX", 1, 2), Op("CX", 2, 3),
]
conditional_measurement = measure(apply(to_rho(zero_state(3)), conditional_ops); label="Conditional circuit")
save_plot("measure_mid.png", conditional_measurement)
Random.seed!(108)
random_ops_example = random_clifford(5, 20; measure_prob=0.2, measure_basis=["MX", "MZ"])
save_plot("random_circuit.png", random_ops_example; fold=12)
println("Basic example figures regenerated.")

## Example 9: local MPS shadows compared with an exact reduced state.
shadow_ops = [Op("H", 1), Op("CX", 1, 2), Op("CX", 2, 3)]
kept_qubits = [1, 3]
# Only this small reference calculation uses a statevector, never the shadows.
rho_exact = partial_trace(apply(zero_state(3), shadow_ops), kept_qubits)
@assert rho_exact ≈ Diagonal([0.5, 0.0, 0.0, 0.5])
shot_counts = [100, 300, 1000, 3000]
repetitions = 3
shadow_errors = zeros(length(shot_counts), repetitions)
Random.seed!(2026) # Reproducible teaching example; independent runs below.
for (i, shots) in enumerate(shot_counts), repetition in 1:repetitions
    rho_hat = shadow(shadow_ops, shots; qubits=kept_qubits,
                     backend=:mps, cutoff=0.0, maxdim=2)
    shadow_errors[i, repetition] = norm(rho_hat - rho_exact) # Frobenius norm.
end
rms_error = sqrt.(vec(mean(shadow_errors.^2; dims=2)))
expected_rms = sqrt.((25 - real(tr(rho_exact * rho_exact))) ./ shot_counts)
println("MPS shadow comparison (shots per estimate, RMS over $repetitions runs, expected RMS):")
for i in eachindex(shot_counts)
    println((shots=shot_counts[i], rms_error=rms_error[i], expected_rms=expected_rms[i]))
end
fig, ax = subplots(figsize=(6.4, 3.9))
for repetition in 1:repetitions
    ax.loglog(shot_counts, shadow_errors[:, repetition], "o";
              color="#aab3bf", markersize=4, alpha=0.7,
              label=repetition == 1 ? "Individual estimates" : "_nolegend_")
end
ax.loglog(shot_counts, rms_error, "o-"; color="#3267a8", label="RMS over 3 runs")
ax.loglog(shot_counts, expected_rms, "--"; color="#d18032", label="Expected RMS")
ax.set_xticks(shot_counts)
ax.set_xticklabels(string.(shot_counts))
ax.set_xlabel("Shots per estimate")
ax.set_ylabel("Frobenius error of the two-qubit density matrix")
ax.set_title("Local MPS shadows of a three-qubit GHZ state", loc="left", fontsize=12)
ax.legend(frameon=false)
save_scientific_plot("shadow_mps_convergence.png", fig, ax)

## Example 11: product-formula evolution and its exact small-system reference.
N = 4
terms = [-1.0, "X,X", -0.8, "Z"]
time = 0.6
H = hamiltonian(N, terms)
exact_state = evolve(zero_state(N), H, time)
trotter_ops = hamiltonian_exp(N, time, terms; dt=0.05)
trotter_state = apply(zero_state(N), trotter_ops)
trotter_measurement = measure(trotter_state; label="Trotter dt=0.05")
exact_measurement = measure(exact_state; label="Exact")
finer_state = apply(zero_state(N), hamiltonian_exp(N, time, terms; dt=0.025))
@assert norm(finer_state-exact_state) < norm(trotter_state-exact_state)
save_plot("trotter_circuit.png", get_layers(trotter_ops)[1:16]; fold=6)
save_plot("trotter_comparison.png", [trotter_measurement, exact_measurement])
fig, ax = subplots(figsize=(5.4, 3.3))
ax.plot(1:N, trotter_measurement.expect, "o-", color="#3267a8", label="Trotter dt=0.05")
ax.plot(1:N, exact_measurement.expect, "x--", color="#d18032", label="Exact")
ax.set_xticks(1:N)
ax.set_xlabel("Qubit")
ax.set_ylabel("Local Z expectation")
ax.set_title("Ising evolution at time 0.6", loc="left", fontsize=12)
ax.legend(frameon=false)
save_scientific_plot("trotter_expect.png", fig, ax)

## Example 12: one finite-shot twirling comparison, with fresh frames per shot.
Random.seed!(112)
rotation_error = 0.2
noise = NoiseModel(false, Noise2("rot_z", rotation_error))
twirl_ops = [Op("H", 1), Op("H", 2), Op("CX", 1, 2),
             Op("CX", 2, 1), Op("H", 1), Op("H", 2)]
twirl_options = [
    Options(circuit_name="Ideal"),
    Options(circuit_name="Coherent noise", noise=noise),
    Options(circuit_name="Randomized frames", noise=noise, twirl=true),
]
twirl_measurements = [measure(twirl_ops, 4000; options=options)
                      for options in twirl_options]
save_plot("circuit_noisy.png", twirl_ops)
save_plot("twirl_measure.png", twirl_measurements)

## Example 14: readout inversion with a known independent bit-flip response.
Random.seed!(114)
readout_ops = [Op("H", 1), Op("CX", 1, 2)]
readout_error = Noise1("bit_flip", 0.1)
readout_options = [
    Options(circuit_name="Ideal"),
    Options(circuit_name="Readout error", readout_noise=readout_error),
    Options(circuit_name="Inverse corrected",
        readout_noise=readout_error, measurement_mitigate=true),
]
readout_measurements = [measure(readout_ops, 4000; options=options, backend=:density_matrix)
                       for options in readout_options]
save_plot("measurement_mitigate.png", readout_measurements)

## Examples 15–17: occupancy-grid routing, Hadamard test, and controlled swap.
layout = Layout(sa.sparse([1 1 1; 1 0 1]))
nonlocal_ops = [Op("H", 1), Op("CX", 1, 5)]
local_ops = [Op("H", 1), Op("CX", 1, 3), Op("I", 5)]
controlled_ops = [Op("H", 1), Op("X", 5; control=1)]
routed_nonlocal = route(nonlocal_ops, layout)
routed_local = route(local_ops, layout)
routed_controlled = route(controlled_ops, layout)
@assert apply(zero_state(5), routed_nonlocal) ≈ apply(zero_state(5), nonlocal_ops)
@assert apply(zero_state(5), routed_controlled) ≈ apply(zero_state(5), controlled_ops)
save_plot("layout1.png", routed_nonlocal)
save_plot("layout2.png", routed_local)
save_plot("layout3.png", routed_controlled)
hadamard_ops = [Op("X", 2), Op("H", 1), Op("Z", 2; control=1), Op("H", 1)]
@assert isapprox(expect(apply(zero_state(2), hadamard_ops), "Z")[1], -1; atol=1e-12)
save_plot("hadamard_test.png", hadamard_ops)
swap_ops = [Op("H", 3), Op("H", 1), Op("SWAP", 2, 3; control=1), Op("H", 1)]
@assert isapprox(expect(apply(zero_state(3), swap_ops), "Z")[1], 0.5; atol=1e-12)
save_plot("swap.png", swap_ops)

## QML: a short variational optimization with an exact energy reference.
N = 3
H = hamiltonian(N, [-1.0, "Z,Z", -0.9, "X"])
exact_energy = eigmin(Hermitian(Matrix(H)))
ops = QuantumOps[Op("RY", _RY, q) for q in 1:N]
append!(ops, [Op("CX", q, q+1) for q in 1:N-1])
energy(parameters) = real(expect(variational_apply(parameters, ops), H))
initial = fill(0.3, nparameters(ops))
result = train(energy, initial; iterations=40)
@assert result.loss >= exact_energy - 1e-10
@assert result.loss <= energy(initial) + 1e-10
save_plot("vqe_c.png", ops)
fig, ax = subplots(figsize=(5.4, 3.3))
ax.plot(1:length(result.history), result.history, "o-", markersize=3,
        color="#3267a8", label="Variational circuit")
ax.axhline(exact_energy, color="#d18032", linestyle="--", label="Exact ground state")
ax.set_xlabel("History entry")
ax.set_ylabel("Energy")
ax.set_title("Three-site Ising model", loc="left", fontsize=12)
ax.legend(frameon=false)
save_scientific_plot("vqe_ising.png", fig, ax)
custom_ops = QuantumOps[Op("local RY", theta -> _RY(theta), 1),
                       Op("CX", 1, 2), Op("local RZ", theta -> _RZ(theta), 2)]
save_plot("vqe_c2.png", custom_ops)
println("Hamiltonian, noise, routing and VQE figures regenerated.")

## QML tutorial: phase regression with a separate held-out grid.
regression_ops = [Op("RY", _RY, 1)]
xtrain = collect(range(-pi, pi; length=13))[1:12]
ytrain = cos.(xtrain .+ 0.65)
xtest = xtrain .+ pi/12
ytest = cos.(xtest .+ 0.65)
function predict(parameters, x)
    encoded = _RY(x) * zero_state(1)
    return real(expect(variational_apply(parameters, regression_ops; init=encoded), gate.Z))
end
loss(parameters) = sum((predict(parameters, x)-y)^2
                       for (x,y) in zip(xtrain,ytrain)) / length(xtrain)
fit = train(loss, [0.2]; iterations=30)
predictions = [predict(fit.parameters, x) for x in xtest]
test_mse = sum(abs2, predictions-ytest) / length(ytest)
@assert test_mse < 1e-10
fig, axes = subplots(1, 2; figsize=(9.0, 4.2))
axes[1].semilogy(1:length(fit.history), max.(fit.history, eps(Float64)), "o-",
                 color="#3267a8", markersize=4)
axes[1].set_xlabel("History entry")
axes[1].set_ylabel("Training mean squared error")
axes[1].set_title("Fit one phase", loc="left", fontsize=12)
xgrid = collect(range(-pi, pi; length=200))
axes[2].plot(xgrid, cos.(xgrid .+ 0.65), color="#8b97a8", linewidth=3,
             alpha=0.55, label="Reference")
axes[2].plot(xgrid, [predict(fit.parameters,x) for x in xgrid], color="#3267a8",
             linewidth=1.5, label="Trained model")
axes[2].scatter(xtrain, ytrain, marker="o", facecolors="none", edgecolors="#3267a8",
                s=35, label="Training labels")
axes[2].scatter(xtest, predictions, marker="x", color="#d18032", s=32,
                label="Held-out predictions")
axes[2].set_xlabel("Encoded angle x (radians)")
axes[2].set_ylabel("Z expectation")
axes[2].set_title("Test on unseen inputs", loc="left", fontsize=12)
axes[2].legend(frameon=false, fontsize=8, loc="upper center",
                bbox_to_anchor=(0.5, -0.2), ncol=2)
for ax in axes
    for side in ("top", "right")
        ax.spines[side].set_visible(false)
    end
    ax.grid(axis="y", color="#e5e9ef", linewidth=0.8)
    ax.set_axisbelow(true)
end
fig.tight_layout(pad=1.1)
savefigure(fig, joinpath(figure_directory, "qml_training.png"); dpi=150)
close(fig)
println("QML held-out mean squared error: $test_mse")

## QML: supervised two-feature XOR classification.

# XOR labels: opposite signs are class 1; equal signs are class 0.
train_levels = [-0.8, -0.3, 0.3, 0.8]
test_levels = collect(range(-0.95, 0.95; length=20))
train_x = [[a,b] for a in train_levels for b in train_levels]
test_x = [[a,b] for a in test_levels for b in test_levels]
label(x) = Int(xor(x[1] > 0, x[2] > 0))
train_y, test_y = label.(train_x), label.(test_x)
@assert isempty(intersect(train_levels, test_levels))

# Encode each feature as an RY angle; precompute these fixed input states.
encoding = [Op("RY", _RY, 1), Op("RY", _RY, 2)]
encode(x) = variational_apply((pi/2) .* x, encoding)
train_states, test_states = encode.(train_x), encode.(test_x)
ansatz = [Op("RY", _RY, 1), Op("RY", _RY, 2),
          Op("CX", 1, 2), Op("RY", _RY, 2)]
function probability(parameters, state)
    psi = variational_apply(parameters, ansatz; init=state)
    return abs2(psi[2]) + abs2(psi[4]) # P(q2=1), basis 00,01,10,11
end
function loss(parameters)
    mean(eachindex(train_states)) do i
        p = clamp(probability(parameters, train_states[i]), 1e-12, 1-1e-12)
        y = train_y[i]
        -y*log(p) - (1-y)*log1p(-p)
    end
end
initial = [0.2, 0.4, -0.1]
result = train(loss, initial; optimizer=:lbfgs, iterations=80)
train_probabilities = probability.(Ref(result.parameters), train_states)
test_probabilities = probability.(Ref(result.parameters), test_states)
train_accuracy = mean((train_probabilities .>= 0.5) .== train_y)
test_accuracy = mean((test_probabilities .>= 0.5) .== test_y)

# An independent attainable solution checks the encoding and readout conventions.
reference(x) = (1 - sin(pi*x[1]/2)*sin(pi*x[2]/2))/2
analytic_parameters = [pi/2, pi/2, 0.0]
@assert all(isapprox(probability(analytic_parameters, encode(x)), reference(x);
                     atol=1e-12) for x in test_x)
@assert test_accuracy >= 0.98
println("Train/test examples: ", (length(train_x), length(test_x)))
println("Initial/final cross-entropy: ", (loss(initial), result.loss))
println("Train/test accuracy: ", (train_accuracy, test_accuracy))
println("Parameters: ", result.parameters)
println("Accepted loss values: ", length(result.history))
println("Held-out cross-entropy: ", mean(-y*log(p)-(1-y)*log1p(-p)
    for (p,y) in zip(test_probabilities,test_y)))

feature_grid = collect(range(-1,1; length=81))
probability_grid = [probability(result.parameters,encode([x,y]))
                    for y in feature_grid, x in feature_grid]
fig, ax = subplots(figsize=(6.6,5.6), dpi=140)
heatmap = ax.imshow(probability_grid; origin="lower", extent=(-1,1,-1,1),
                    vmin=0, vmax=1, cmap="coolwarm", interpolation="bilinear")
# The fitted angles match the analytic reference; its exact boundary is x₁x₂=0.
ax.axvline(0; color="#263346", linestyle="--", linewidth=1.4)
ax.axhline(0; color="#263346", linestyle="--", linewidth=1.4)
for (class, marker, face) in ((0,"o","white"),(1,"^","#263346"))
    points = train_x[train_y .== class]
    ax.scatter(first.(points),last.(points); marker=marker, s=58,
               facecolors=face, edgecolors="#263346", linewidths=1.1,
               label="Training class $class", zorder=3)
end
ax.plot([],[]; color="#263346", linestyle="--", linewidth=1.4, label="Reference p = 0.5")
ax.set_xlabel("Feature x₁", fontsize=11)
ax.set_ylabel("Feature x₂", fontsize=11)
ax.set_xticks([-1,-0.5,0,0.5,1])
ax.set_yticks([-1,-0.5,0,0.5,1])
ax.set_title("Learned XOR probability", loc="left", fontsize=13, pad=12)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.16),ncol=3,
          frameon=false,fontsize=9,handlelength=1.7,columnspacing=1.1)
colorbar = fig.colorbar(heatmap; ax=ax, fraction=0.046, pad=0.04)
colorbar.set_label("P(q₂ = 1): class 1 probability", fontsize=10, labelpad=10)
fig.tight_layout(pad=1.1)
savefigure(fig,joinpath(figure_directory, "qml_classifier.png");dpi=160)
close(fig)

## Example 19: finite-size monitored trajectories, without deliberate truncation.
Random.seed!(119)
unitary_ops = random_ops(6, 4)
dense = apply(zero_state(6), unitary_ops)
psi = apply(zero_state(N_MPS(6)), unitary_ops; cutoff=0.0, maxdim=8)
@assert isapprox(abs2(dot(dense, to_state(psi))), 1; atol=1e-10)
illustrative_ops = random_ops(6, 3; measure_prob=0.2, measure_basis=["MZ"])
save_plot("random_ops_m.png", illustrative_ops; fold=12)

sizes = [4, 6]
probabilities = [0.0, 0.25, 0.5, 0.75, 1.0]
depth = 8
trajectories = 24
entropy_samples = Dict{Tuple{Int,Float64},Vector{Float64}}()
for N in sizes, p in probabilities
    samples = Float64[]
    for repetition in 1:trajectories
        ops = random_ops(N, depth; measure_prob=p, measure_basis=["MZ"])
        trajectory_state = apply(zero_state(N_MPS(N)), ops; cutoff=0.0, maxdim=2^(N÷2))
        push!(samples, entanglement_entropy(trajectory_state, N÷2))
    end
    @assert all(S -> -1e-10 <= S <= (N÷2)*log(2) + 1e-10, samples)
    entropy_samples[N, p] = samples
    println("N=$N, p=$p: mean entropy=$(round(mean(samples); digits=4)) ± $(round(std(samples)/sqrt(trajectories); digits=4)) nats (SE)")
end
@assert all(N -> maximum(abs, entropy_samples[N, 1.0]) < 1e-10, sizes)
fig, ax = subplots(figsize=(6.1, 3.9))
for (i, N) in enumerate(sizes)
    means = [mean(entropy_samples[N,p]) for p in probabilities]
    errors = [std(entropy_samples[N,p])/sqrt(trajectories) for p in probabilities]
    ax.errorbar(probabilities, means; yerr=errors, marker=i==1 ? "o" : "s",
                color=i==1 ? "#3267a8" : "#d18032", linewidth=1.3, capsize=3, label="N=$N")
end
ax.set_xlabel("Z-measurement probability per site per layer")
ax.set_ylabel("Mean half-chain trajectory entropy (nats)")
ax.set_title("8 layers · 24 trajectories per point", loc="left", fontsize=12)
ax.set_xticks(probabilities)
ax.legend(frameon=false)
save_scientific_plot("mie.png", fig, ax)
println("All documentation figures regenerated.")
