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
    for storage in (identity, sa.sparse)
        state = storage(psi)
        for q in 1:4
            reference = projected_partial_trace(rho, fill(2, 4), setdiff(1:4, [q]))
            @test partial_trace(state, q) ≈ reference
        end
        for q in 1:3
            reference = projected_partial_trace(rho, fill(2, 4), setdiff(1:4, [q, q + 1]))
            @test partial_trace(state, q, q + 1) ≈ reference
            @test partial_trace(state, q + 1, q) ≈ reference
        end
        @test BlueTangle.bipartition_trace(storage(rho)) ≈
              projected_partial_trace(rho, fill(2, 4), [1, 2])

        for keep in ([0], [5], [1, 5], [1, 1], [1.0])
            @test_throws ArgumentError partial_trace(state, keep)
        end
        for q in (0, 5)
            @test_throws ArgumentError partial_trace(state, q)
        end
        for q in 1:4
            @test_throws ArgumentError partial_trace(state, q, q)
        end
        @test_throws ArgumentError partial_trace(state, 1, 3)
        @test_throws ArgumentError partial_trace(state, 0, 1)
        @test_throws ArgumentError partial_trace(state, 4, 5)
    end
    for state in (ComplexF64[], ones(ComplexF64, 3))
        @test_throws DimensionMismatch partial_trace(state, Int[])
    end

    # Known product states expose bit-order and conjugation errors directly.
    plus_y = ComplexF64[1, im] / sqrt(2)
    product = kron(ComplexF64[0, 1], plus_y, ComplexF64[1, 0])
    @test partial_trace(product, 1) ≈ [0 0; 0 1]
    @test partial_trace(product, 2) ≈ plus_y * plus_y'
    @test partial_trace(product, 3) ≈ [1 0; 0 0]
    @test partial_trace(2product, Int[]) ≈ fill(4, 1, 1)
    @test partial_trace(plus_y, 1) ≈ plus_y * plus_y'
    @test partial_trace(ComplexF64[2im], Int[]) == fill(4, 1, 1)

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
    for storage in (identity, sa.sparse)
        @test BlueTangle.partial_trace_second_subsystem(storage(mixed), 2, 6) ≈
              projected_partial_trace(mixed, [2, 6], [2])
        @test_throws DimensionMismatch BlueTangle.partial_trace_second_subsystem(storage(mixed), 2, 2)
        @test_throws DimensionMismatch BlueTangle.bipartition_trace(storage(ones(4, 3)))
        @test_throws ArgumentError BlueTangle.bipartition_trace(storage(ones(8, 8)))
    end
    for invalid_dims in ([0, 6], [-2, 6], [2.0, 6.0])
        @test_throws ArgumentError partial_trace(mixed, invalid_dims, [1])
    end
end
