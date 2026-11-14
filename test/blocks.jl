@testset "Static blocks preserve operation semantics" begin
    damping = NoiseModel(Noise1("amplitude_damping", 1.0), false)
    states = (zero_state(1), zero_state(N_MPS(1)), ComplexF64[1 0; 0 0])
    comparable(state) = state isa it.MPS ? to_state(state) : state

    for initial in states
        operations = QuantumOps[Op("X", 1)]
        block = Op("X block", operations)
        expected = apply(operations, initial; noise=damping)
        @test comparable(apply(initial, block; noise=damping)) ≈ comparable(expected)
        @test comparable(variational_apply(Float64[], QuantumOps[block]; init=initial, noise=damping)) ≈ comparable(expected)
        @test comparable(expected) ≈ comparable(initial)

        protected_operations = QuantumOps[Op("X", 1; noisy=false)]
        protected = Op("protected X", protected_operations)
        protected_expected = apply(protected_operations, initial; noise=damping)
        @test comparable(apply(initial, protected; noise=damping)) ≈ comparable(protected_expected)
        @test comparable(variational_apply(Float64[], QuantumOps[protected]; init=initial, noise=damping)) ≈ comparable(protected_expected)
        @test expect(protected_expected, Op("Z", 1)) ≈ -1.0

        nested = Op("nested X block", QuantumOps[block])
        @test comparable(apply(initial, nested; noise=damping)) ≈ comparable(expected)
        @test comparable(variational_apply(Float64[], QuantumOps[nested]; init=initial, noise=damping)) ≈ comparable(expected)

        recorded_operations = QuantumOps[Op("X", 1; noisy=false), Op("MZ", 1)]
        measured = Op("measured block", recorded_operations)
        nested_measured = Op("nested measured block", QuantumOps[measured])
        expected_state, expected_bits = apply(recorded_operations, initial; noise=damping, track_measurements=true)
        @test expected_bits == [1]
        for measured_block in (measured, nested_measured)
            result_state, result_bits = apply(initial, measured_block; noise=damping, track_measurements=true)
            @test comparable(result_state) ≈ comparable(expected_state)
            @test result_bits == expected_bits
        end

        # Opaque matrix/callback blocks retain their existing behavior: there
        # is no constituent gate list on which to attach local gate noise.
        for opaque in (Op("matrix block", gate.X), Op("callback block", state -> apply(state, Op("X", 1))))
            @test expect(apply(initial, opaque; noise=damping), Op("Z", 1)) ≈ -1.0
        end
    end
end

@testset "Captured parameters inside static blocks" begin
    theta = [1.1]
    gamma = 0.4
    noise = NoiseModel(Noise1("amplitude_damping", gamma), false)
    rho0 = ComplexF64[1 0; 0 0]
    for depth in 1:2
        function block_loss(p; init=zero_state(1), noise=false)
            operations = QuantumOps[Op("matrix", _RY(p[1]), 1)]
            for _ in 1:depth
                operations = QuantumOps[Op("static block", operations)]
            end
            return expect(variational_apply(Float64[], operations; init=init, noise=noise), gate.Z)
        end
        # Materialized gates still carry AD parameters through fixed blocks.
        @test only(BlueTangle.ForwardDiff.gradient(block_loss, theta)) ≈ -sin(only(theta))
        exact_loss(p) = block_loss(p; init=rho0, noise=noise)
        @test exact_loss(theta) ≈ gamma + (1-gamma)*cos(only(theta))
        @test only(BlueTangle.ForwardDiff.gradient(exact_loss, theta)) ≈ -(1-gamma)*sin(only(theta))
        @test_throws ArgumentError BlueTangle.ForwardDiff.gradient(p -> block_loss(p; noise=noise), theta)
        @test_throws ArgumentError BlueTangle.ForwardDiff.gradient(p -> block_loss(p; init=zero_state(N_MPS(1))), theta)
    end
    function conditional_loss(p)
        conditional = ifOp("MZ", 1, Op("matrix", _RY(p[1]), 1), Op("I", 1))
        return expect(variational_apply(Float64[], QuantumOps[conditional]; N=1), gate.Z)
    end
    @test_throws ArgumentError BlueTangle.ForwardDiff.gradient(conditional_loss, theta)
end

@testset "MPS block truncation and parameter boundaries" begin
    initial = zero_state(N_MPS(2))
    operations = QuantumOps[Op("RY(1.1)", 1), Op("RY(0.8)", 2), Op("CX", 1, 2)]
    block = Op("entangling block", operations)
    nested = Op("nested entangling block", QuantumOps[block])
    expected = apply(operations, initial; cutoff=0.5, maxdim=1)
    for wrapped in (block, nested)
        result = apply(initial, wrapped; cutoff=0.5, maxdim=1)
        @test to_state(result) ≈ to_state(expected)
        @test maxlinkdim(result) == 1
        variational_result = variational_apply(Float64[], QuantumOps[wrapped];
            init=initial, cutoff=0.5, maxdim=1)
        @test to_state(variational_result) ≈ to_state(expected)
    end
    @test to_state(initial) ≈ zero_state(2)

    # Parameters remain attached to top-level local gates. A named static
    # block does not introduce recursive variational parameter binding.
    parameterized_block = Op("parameterized block", QuantumOps[Op("RY", _RY, 1)])
    @test nparameters(QuantumOps[parameterized_block]) == 0
    @test_throws DimensionMismatch variational_apply([0.3], QuantumOps[parameterized_block])
end
