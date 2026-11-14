@testset "State-independent Kraus sampling for MPS" begin
    logical = normalize(ComplexF64[1, 0.3im, -0.4, 0.2+0.7im, 0.1, 0.5im, -0.2im, 0.3])
    mps = to_MPS(logical, N_MPS(3); cutoff=0.0)
    # A correlated random-unitary channel on nonadjacent, reversed-order wires
    # checks both physical ordering and state-independent branch probabilities.
    kraus = [sqrt(0.37) * kron(gate.X, gate.Z),
             sqrt(0.63) * kron(gate.H, gate.S)]
    channel = custom_noise(2, "correlated unitary", kraus)
    op = Op(channel, 3, 1)
    observed = Set{Int}()
    for seed in 1:24
        Random.seed!(seed)
        expected = apply(logical, op)
        Random.seed!(seed)
        actual = to_state(apply(mps, op))
        @test isapprox(actual, expected; atol=1e-10)
        push!(observed, abs2(dot(expected, normalize(Op("first", kraus[1], 3, 1).expand(3)*logical))) > 1-1e-10 ? 1 : 2)
    end
    @test observed == Set([1,2])

    # Global phase must survive the identity shortcut; an identity with
    # explicit truncation settings must still use the ordinary apply path.
    phase = custom_noise(2, "phase identity", [im*Matrix{ComplexF64}(I, 4, 4)])
    @test to_state(apply(mps, Op(phase, 1, 3))) ≈ im*logical
    identity = custom_noise(2, "identity", [Matrix{Float64}(I, 4, 4)])
    for settings in ((cutoff=0.0,), (maxdim=1,))
        actual = apply(mps, Op(identity, 1, 3); settings...)
        sites = BlueTangle.get_M(mps)
        expected = normalize(it.apply(it.op(gate.I ⊗ gate.I, sites[3], sites[1]), mps; settings...))
        @test to_state(actual) ≈ to_state(expected)
    end

    # State-dependent channels must keep the original Born-probability path.
    for gamma in (0.4, 1e-9)
        damping = OpQC("amplitude_damping", gamma, 2)
        for seed in 1:8
            Random.seed!(seed)
            expected = apply(logical, damping)
            Random.seed!(seed)
            actual = to_state(apply(mps, damping))
            @test isapprox(actual, expected; atol=1e-12)
        end
    end

    # Two full-rank weak filters can both look unitary under default ≈.
    # Place the known RNG draw between their average and actual Born weights:
    # using state-independent probabilities would select the wrong branch.
    Random.seed!(189)
    threshold = rand()
    delta = 1e-9
    polarization = real(dot(logical, Op("Z", 2).expand(3)*logical))
    p = threshold - delta*polarization/2
    filters = [diagm(sqrt.([p+delta, p-delta])),
               diagm(sqrt.([1-p-delta, 1-p+delta]))]
    weak_filter = custom_noise(1, "weak filters", filters)
    Random.seed!(189)
    expected = apply(logical, Op(weak_filter, 2))
    Random.seed!(189)
    actual = to_state(apply(mps, Op(weak_filter, 2)))
    @test isapprox(actual, expected; atol=1e-12, rtol=1e-12)

    # Do not cache classification by channel name or identity: supported direct
    # edits of its shared Kraus list must be reflected at the next application.
    mutable_channel = custom_noise(1, "mutable", [gate.I])
    empty!(mutable_channel.kraus)
    append!(mutable_channel.kraus, noise_model1("amplitude_damping", 0.4))
    for seed in 1:8
        Random.seed!(seed)
        expected = apply(logical, Op(mutable_channel, 2))
        Random.seed!(seed)
        actual = to_state(apply(mps, Op(mutable_channel, 2)))
        @test isapprox(actual, expected; atol=1e-12)
    end
end
