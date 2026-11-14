@testset "Partial trace" begin
    # An independent reference: project each traced subsystem onto its basis,
    # then sum the resulting density matrices of the retained subsystems.
    function projected_partial_trace(rho, dims, trace_out)
        reduced_dim = prod(dims[setdiff(eachindex(dims), trace_out)])
        reduced = zeros(ComplexF64, reduced_dim, reduced_dim)
        for basis in Iterators.product((1:dims[q] for q in trace_out)...)
            factors = Matrix{ComplexF64}[]
            for q in eachindex(dims)
                traced_position = findfirst(==(q), trace_out)
                if isnothing(traced_position)
                    push!(factors, Matrix{ComplexF64}(I, dims[q], dims[q]))
                else
                    bra = zeros(ComplexF64, 1, dims[q])
                    bra[basis[traced_position]] = 1
                    push!(factors, bra)
                end
            end
            projector = reduce(kron, factors)
            reduced += projector * rho * projector'
        end
        return reduced
    end

    rng = MersenneTwister(42)
    psi = normalize(randn(rng, ComplexF64, 16))
    rho = psi * psi'
    for mask in 0:15
        keep = [q for q in 1:4 if !iszero(mask & (1 << (q - 1)))]
        trace_out = setdiff(1:4, keep)
        reference = projected_partial_trace(rho, fill(2, 4), trace_out)
        @test partial_trace(psi, keep) ≈ reference
        @test partial_trace(sa.sparse(psi), keep) ≈ reference
        @test partial_trace(rho, fill(2, 4), trace_out) ≈ reference
        @test partial_trace(sa.sparse(rho), fill(2, 4), trace_out) ≈ reference
    end
    @test partial_trace(psi, [4, 1]) ≈ partial_trace(psi, [1, 4])
    for q in 1:4
        @test partial_trace(psi, q) ≈ partial_trace(psi, [q])
    end
    for q in 1:3
        @test partial_trace(psi, q, q + 1) ≈ partial_trace(psi, [q, q + 1])
    end

    # Mixed states and unequal subsystem dimensions exercise the matrix API.
    dims = [2, 3, 2]
    amplitudes = randn(rng, ComplexF64, 12, 5)
    mixed = amplitudes * amplitudes'
    mixed /= tr(mixed)
    for trace_out in (Int[], [1], [2], [1, 3], [1, 2, 3])
        reference = projected_partial_trace(mixed, dims, trace_out)
        @test partial_trace(mixed, dims, trace_out) ≈ reference
        @test partial_trace(sa.sparse(mixed), dims, trace_out) ≈ reference
    end
    @test partial_trace(mixed, dims, [3, 1]) ≈ partial_trace(mixed, dims, [1, 3])
    @test_throws DimensionMismatch partial_trace(mixed, [2, 2], [1])
    @test_throws ArgumentError partial_trace(mixed, dims, [1, 1])
    @test_throws ArgumentError partial_trace(mixed, dims, [4])
end
