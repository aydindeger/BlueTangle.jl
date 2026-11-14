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
ideal = measure(apply(bell_ops, rho0); label="Ideal")
noisy = measure(apply(bell_ops, rho0; noise=NoiseModel("depolarizing", 0.08)); label="With gate noise")
save_plot("measurement_comparison.png", [ideal, noisy], ["Ideal", "With gate noise"])

## Examples 6 and 8: conditional measurement and a random Clifford circuit.
conditional_ops = QuantumOps[
    Op("X", 1), Op("MX", 1; if0="Z", if1="H"),
    Op("CX", 1, 2), Op("CX", 2, 3),
]
conditional_measurement = measure(to_rho(compile(conditional_ops)); label="Conditional circuit")
save_plot("measure_mid.png", conditional_measurement)
Random.seed!(108)
random_circuit = compile(random_clifford(5, 20; measure_prob=0.2, measure_basis=["MX", "MZ"]))
save_plot("random_circuit.png", random_circuit; fold=12)
println("Basic example figures regenerated.")

## Example 11: product-formula evolution and its exact small-system reference.
N = 4
terms = [-1.0, "X,X", -0.8, "Z"]
time = 0.6
H = hamiltonian(N, terms)
exact_state = evolve(zero_state(N), H, time)
trotter_ops = hamiltonian_exp(N, time, terms; dt=0.05)
trotter_circuit = compile(trotter_ops)
trotter_state = to_state(trotter_circuit)
trotter_measurement = measure(trotter_state; label="Trotter dt=0.05")
exact_measurement = measure(exact_state; label="Exact")
finer_state = to_state(compile(hamiltonian_exp(N, time, terms; dt=0.025)))
@assert norm(finer_state-exact_state) < norm(trotter_state-exact_state)
save_plot("trotter_circuit.png", trotter_circuit.layers[1:16]; fold=6)
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
ideal_circuit = compile(twirl_ops, Options(circuit_name="Ideal"))
noisy_circuit = compile(twirl_ops, Options(circuit_name="Coherent noise", noise=noise))
twirled_circuit = compile(twirl_ops,
    Options(circuit_name="Randomized frames", noise=noise, twirl=true))
twirl_measurements = [measure(circuit, 4000) for circuit in
                      (ideal_circuit, noisy_circuit, twirled_circuit)]
save_plot("circuit_noisy.png", noisy_circuit)
save_plot("twirl_measure.png", twirl_measurements)

## Example 14: readout inversion with a known independent bit-flip response.
Random.seed!(114)
readout_ops = [Op("H", 1), Op("CX", 1, 2)]
readout_error = Noise1("bit_flip", 0.1)
readout_circuits = [
    compile(readout_ops, Options(circuit_name="Ideal")),
    compile(readout_ops, Options(circuit_name="Readout error", readout_noise=readout_error)),
    compile(readout_ops, Options(circuit_name="Inverse corrected",
        readout_noise=readout_error, measurement_mitigate=true)),
]
readout_measurements = [measure(circuit, 4000; backend=:density_matrix)
                       for circuit in readout_circuits]
save_plot("measurement_mitigate.png", readout_measurements)

## Examples 15–17: occupancy-grid routing, Hadamard test, and controlled swap.
layout = Layout(sa.sparse([1 1 1; 1 0 1]))
nonlocal_ops = [Op("H", 1), Op("CX", 1, 5)]
local_ops = [Op("H", 1), Op("CX", 1, 3), Op("I", 5)]
controlled_ops = [Op("H", 1), Op("X", 5; control=1)]
nonlocal_circuit = compile(nonlocal_ops; layout=layout)
local_circuit = compile(local_ops; layout=layout)
controlled_circuit = compile(controlled_ops; layout=layout)
@assert to_state(nonlocal_circuit) ≈ apply(nonlocal_ops, zero_state(5))
@assert to_state(controlled_circuit) ≈ apply(controlled_ops, zero_state(5))
save_plot("layout1.png", nonlocal_circuit)
save_plot("layout2.png", local_circuit)
save_plot("layout3.png", controlled_circuit)
hadamard_ops = [Op("X", 2), Op("H", 1), Op("Z", 2; control=1), Op("H", 1)]
@assert isapprox(expect(apply(hadamard_ops, zero_state(2)), "Z")[1], -1; atol=1e-12)
save_plot("hadamard_test.png", hadamard_ops)
swap_ops = [Op("H", 3), Op("H", 1), Op("SWAP", 2, 3; control=1), Op("H", 1)]
@assert isapprox(expect(apply(swap_ops, zero_state(3)), "Z")[1], 0.5; atol=1e-12)
save_plot("swap.png", swap_ops)

## Example 18: a short variational optimization with an exact energy reference.
N = 3
H = hamiltonian(N, [-1.0, "Z,Z", -0.9, "X"])
exact_energy = eigmin(Hermitian(Matrix(H)))
ansatz = AnsatzOptions(N=N, ops=["RY", "CX"], loss=H,
    deep_circuit=true, model="lbfgs", paramshift=true,
    pars_initial=fill(0.3, 100), number_of_iterations=40)
energies, parameters = VQE(ansatz)
@assert real(expect(variational_apply(parameters, ansatz), H)) >= exact_energy - 1e-10
save_plot("vqe_c.png", ansatz)
fig, ax = subplots(figsize=(5.4, 3.3))
ax.plot(1:length(energies), energies, "o-", markersize=3, color="#3267a8", label="Variational circuit")
ax.axhline(exact_energy, color="#d18032", linestyle="--", label="Exact ground state")
ax.set_xlabel("Objective evaluation")
ax.set_ylabel("Energy")
ax.set_title("Three-site Ising model", loc="left", fontsize=12)
ax.legend(frameon=false)
save_scientific_plot("vqe_ising.png", fig, ax)
custom_ops = QuantumOps[Op("local RY", theta -> _RY(theta), 1),
                       Op("CX", 1, 2), Op("local RZ", theta -> _RZ(theta), 2)]
custom_ansatz = AnsatzOptions(N=2, ops=custom_ops,
    loss=hamiltonian(2, [-1.0, "Z,Z", -0.7, "X"]), pars_initial=[0.2, 0.4])
save_plot("vqe_c2.png", custom_ansatz)
println("Hamiltonian, noise, routing and VQE figures regenerated.")

## Example 19: finite-size monitored trajectories, without deliberate truncation.
Random.seed!(119)
unitary_ops = random_ops(6, 4)
dense = apply(unitary_ops, zero_state(6))
psi = apply(unitary_ops, zero_state(N_MPS(6)); cutoff=0.0, maxdim=8)
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
        trajectory_state = apply(ops, zero_state(N_MPS(N)); cutoff=0.0, maxdim=2^(N÷2))
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
println("All 20 currently referenced documentation figures regenerated.")
