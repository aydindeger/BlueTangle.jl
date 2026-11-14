@testset "Hamiltonian evolution" begin
    for terms in ([0.7,"X"],[0.7,"Y"],[0.7,"Z"],[0.7,"I"],[0.7,"I,Z"],[0.7,"X,I,Y"])
        N = length(split(terms[2],','))
        H = hamiltonian(N,terms)
        for time in (0.13,-0.13,0.01)
            ops = hamiltonian_exp(N,time,terms;dt=0.05)
            @test sa.sparse(ops) ≈ exp(-im*time*Matrix(H))
        end
    end
    terms = [0.4,"X",0.7,"Z,Z"]
    N = 3
    H = hamiltonian(N,terms)
    expected = evolve(zero_state(N),H,0.23)
    coarse = apply(zero_state(N), hamiltonian_exp(N,0.23,terms;dt=0.05))
    fine = apply(zero_state(N), hamiltonian_exp(N,0.23,terms;dt=0.005))
    @test norm(fine-expected) < norm(coarse-expected)/5
    @test evolve(zero_state(N)*zero_state(N)',H,0.23) ≈ expected*expected'
    @test isempty(hamiltonian_exp(2,0.0,terms))
    for terms in ([0.4,"Y"],[0.4,"Z,Z"],[0.4,"Z,Z,Z"])
        H = hamiltonian((3,2),terms)
        @test sa.sparse(hamiltonian_exp((3,2),0.07,terms;dt=0.03)) ≈ exp(-im*0.07*Matrix(H))
    end
    @test_throws ArgumentError hamiltonian_exp(2,0.1,[0.4,"Q"])
    @test_throws ArgumentError hamiltonian_exp(2,0.1,terms;dt=0)
    @test_throws ArgumentError hamiltonian(2,terms,"invalid")
end

@testset "Product-formula orders" begin
    # Expand with the requested register size: identity padding need not leave
    # a gate on the last qubit, so circuit support alone cannot determine N.
    function circuit_matrix(ops, N)
        U = Matrix{ComplexF64}(I, 2^N, 2^N)
        for op in ops
            U = op.expand(N) * U
        end
        return U
    end

    @testset "Exact Pauli exponentials and time intervals" begin
        for order in (1, 2, 4), word in ("X", "Y", "Z", "I", "I,I", "I,Y", "X,I", "X,I,Y")
            N = length(split(word, ','))
            terms = [0.7, word]
            H = Matrix(hamiltonian(N, terms))
            for time in (0.13, -0.13, 0.01)
                ops = hamiltonian_exp(N, time, terms; dt=0.05, order=order)
                @test circuit_matrix(ops, N) ≈ exp(-im * time * H)
            end
        end

        for order in (1, 2, 4)
            terms = [0.7, "Y"]
            for time in (0.001, -0.001, 0.0)
                direction = time < 0 ? -1 : 1
                ops = hamiltonian_exp(1, time, terms; dt=0.05, order=order, full=false)
                @test circuit_matrix(ops, 1) ≈ exp(-im * direction * 0.05 * 0.7 * gate.Y)
            end
            @test isempty(hamiltonian_exp(2, 0.0, [0.7, "X"]; order=order))
            @test isempty(hamiltonian_exp(2, 0.13, [0.0, "X"]; order=order))
            @test isempty(hamiltonian_exp(2, 0.13, []; order=order))
            @test isempty(hamiltonian_exp(2, 0.13, [0.7, "X,Y,Z"]; order=order))

            # Comparing the off-diagonal entry avoids hiding a lost tiny step
            # behind the order-one diagonal entries in a matrix norm.
            tiny_time = 1e-20
            ops = hamiltonian_exp(1, tiny_time, [0.7, "X"]; dt=0.1, order=order)
            @test !isempty(ops)
            @test isapprox(circuit_matrix(ops, 1)[1, 2], -im * sin(0.7 * tiny_time); rtol=1e-12, atol=0)

            # Large couplings keep the angle visible while a tiny total time
            # exercises the general, noncommuting step scheduler.
            scaled_terms = [1e20, "X", 0.5e20, "Z"]
            full = hamiltonian_exp(1, tiny_time, scaled_terms; dt=0.1, order=order)
            step = hamiltonian_exp(1, tiny_time, scaled_terms; dt=tiny_time, order=order, full=false)
            @test circuit_matrix(full, 1) ≈ circuit_matrix(step, 1)
        end
    end

    @testset "Convergence against dense evolution" begin
        # Test the full unitary, not a special initial state. The translated
        # X,X,Y strings overlap and anticommute, even though they share a label.
        for (N, terms) in ((3, [0.7, "X", -0.4, "Z,Z", 0.2, "Y"]),
                           (4, [0.7, "X,X,Y"]))
            time = 0.4
            exact = exp(-im * time * Matrix(hamiltonian(N, terms)))
            for order in (1, 2, 4)
                coarse = circuit_matrix(hamiltonian_exp(N, time, terms; dt=0.1, order=order), N)
                fine = circuit_matrix(hamiltonian_exp(N, time, terms; dt=0.05, order=order), N)
                @test isapprox(norm(coarse - exact) / norm(fine - exact), 2.0^order; rtol=0.15)
                @test fine' * fine ≈ Matrix{ComplexF64}(I, 2^N, 2^N)
            end
        end

        terms = [0.4, "X", 0.7, "Z,Z"]
        default = circuit_matrix(hamiltonian_exp(3, 0.23, terms; dt=0.05), 3)
        @test default ≈ circuit_matrix(hamiltonian_exp(3, 0.23, terms; dt=0.05, order=1), 3)
        for order in (1, 2, 4), time in (0.23, -0.23)
            exact = exp(-im * time * Matrix(hamiltonian(3, terms)))
            coarse = circuit_matrix(hamiltonian_exp(3, time, terms; dt=0.05, order=order), 3)
            fine = circuit_matrix(hamiltonian_exp(3, time, terms; dt=0.025, order=order), 3)
            @test norm(fine - exact) < norm(coarse - exact)
        end
        for order in (2, 4)
            forward = circuit_matrix(hamiltonian_exp(3, 0.07, terms; dt=0.07, order=order), 3)
            backward = circuit_matrix(hamiltonian_exp(3, -0.07, terms; dt=0.07, order=order), 3)
            @test backward ≈ forward'
        end
        for order in (1, 2, 4)
            ising = BlueTangle.trotter_ising(3, 0.23, 0.7, 0.4; dt=0.05, order=order)
            @test circuit_matrix(ising, 3) ≈ circuit_matrix(hamiltonian_exp(3, 0.23, [-0.4, "X", -0.7, "Z,Z"]; dt=0.05, order=order), 3)
        end
    end

    @testset "Application across simulation backends" begin
        state = normalize(ComplexF64[complex(q, 9-q) for q in 1:8])
        rho = 0.7 * state * state' + 0.3 * Matrix{ComplexF64}(I, 8, 8) / 8
        initial_mps = to_MPS(state, N_MPS(3); cutoff=0.0)
        for order in (1, 2, 4)
            ops = hamiltonian_exp(3, 0.07, [-0.3, "X", 0.2, "Z,Z", 0.1, "Y,I,Z"]; dt=0.05, order=order)
            # Matrix expansion is independent of statevector, density, and
            # tensor application; include a nonlocal Pauli string and complex
            # amplitudes to expose basis-rotation or conjugation mistakes.
            U = circuit_matrix(ops, 3)
            @test apply(state, ops) ≈ U * state
            @test apply(rho, ops) ≈ U * rho * U'
            @test isapprox(to_state(apply(initial_mps, ops; cutoff=0.0)), U * state; atol=1e-10)
        end
    end

    @testset "Two-dimensional lattice enumeration" begin
        # Build the 2×3 Hamiltonian independently from explicitly listed bonds.
        for (style, bonds) in ((:rowwise, [(1,2), (2,3), (4,5), (5,6), (1,4), (2,5), (3,6)]),
                               (:colwise, [(1,3), (3,5), (2,4), (4,6), (1,2), (3,4), (5,6)]))
            H = zeros(ComplexF64, 64, 64)
            for (a, b) in bonds
                factors = [q == a ? gate.X : q == b ? gate.Y : gate.I for q in 1:6]
                H += 0.37 * reduce(kron, factors)
            end
            for a in 1:6
                H -= 0.21 * reduce(kron, [q == a ? gate.Z : gate.I for q in 1:6])
            end
            for (order, tolerance) in ((1, 0.05), (2, 1e-3), (4, 1e-5))
                ops = hamiltonian_exp((2,3), -0.13, [0.37, "X,Y", -0.21, "Z"]; dt=0.04, order=order, enumeration_style=style)
                @test norm(circuit_matrix(ops, 6) - exp(im * 0.13 * H)) < tolerance
            end
        end
    end

    @testset "Gate reductions" begin
        # The commuting X and ZZ blocks share their outer half steps. With
        # three time steps this leaves four X layers and three ZZ layers;
        # fourth order composes nine symmetric stages instead of three.
        for (order, gates, cnots) in ((2, 43, 18), (4, 121, 54))
            ops = hamiltonian_exp(4, 0.3, [0.4, "X", 0.7, "Z,Z"]; dt=0.1, order=order)
            @test length(ops) == gates
            @test count(op -> op.name == "CX", ops) == cnots
        end
        for order in (1, 2, 4)
            ops = hamiltonian_exp(4, 0.3, [0.7, "X,X"]; dt=0.01, order=order)
            @test length(ops) == 21
            @test count(op -> op.name == "CX", ops) == 6
            @test circuit_matrix(ops, 4) ≈ exp(-im * 0.3 * Matrix(hamiltonian(4, [0.7, "X,X"])))
            for word in ("X", "Y", "Z")
                @test length(hamiltonian_exp(1, 0.3, [0.7, word]; dt=0.01, order=order)) == 1
            end
        end
    end

    @testset "Input validation before circuit construction" begin
        for terms in ([0.4, "Q"], [0.4, "Q,Q,Q"], [0.4, "X,"], [0.4, ""], [0.4, 1])
            for time in (0.0, 0.1)
                @test_throws ArgumentError hamiltonian_exp(2, time, terms)
            end
        end
        @test_throws ArgumentError hamiltonian_exp(2, 0.1, [0.4, "X", 0.3])
        for coupling in (NaN, Inf, 0.1im)
            @test_throws ArgumentError hamiltonian_exp(2, 0.0, [coupling, "X"])
        end
        for order in (0, 3, 6)
            @test_throws ArgumentError hamiltonian_exp(2, 0.1, [0.4, "X"]; order=order)
        end
        for dt in (0.0, -0.1, NaN, Inf)
            @test_throws ArgumentError hamiltonian_exp(2, 0.1, [0.4, "X"]; dt=dt)
        end
        for time in (NaN, Inf, -Inf)
            @test_throws ArgumentError hamiltonian_exp(2, time, [0.4, "X"])
        end
        @test_throws ArgumentError hamiltonian_exp((2,3), 0.1, [0.4, "X"]; enumeration_style=:invalid)
    end
end
