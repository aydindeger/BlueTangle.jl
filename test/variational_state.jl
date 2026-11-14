@testset "Local variational statevector kernel" begin
    rng = MersenneTwister(72)
    one_qubit = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
    two_qubit = Matrix(qr(randn(rng, ComplexF64, 4, 4)).Q)
    # Compare arbitrary complex gates with the independently expanded operator,
    # covering both operand orders, nonlocal pairs, and every possible control.
    for N in 1:5
        state = normalize(randn(rng, ComplexF64, 2^N))
        original = copy(state)
        for q in 1:N, control in [-2; setdiff(1:N, [q])]
            op = Op("custom one", one_qubit, q; control=control)
            @test BlueTangle._variational_gate(state, op) ≈ op.expand(N)*state
        end
        for q in 1:N, target in setdiff(1:N, [q]), control in [-2; setdiff(1:N, [q,target])]
            op = Op("custom two", two_qubit, q, target; control=control)
            @test BlueTangle._variational_gate(state, op) ≈ op.expand(N)*state
        end
        @test state == original
    end

    @test BlueTangle._variational_gate(zero_state(3), Op("X",1)) == [0,0,0,0,1,0,0,0]
    @test BlueTangle._variational_gate([1.0,0.0], Op("Y",1)) == [0,im]
    @test BlueTangle._variational_gate(sa.sparse([1.0,0.0]), Op("H",1)) isa Vector

    owned = normalize(randn(rng, ComplexF64, 16))
    fixed = Op("custom two",two_qubit,4,1;control=2)
    expected = fixed.expand(4)*owned
    @test BlueTangle._variational_gate!(owned,fixed) === owned
    @test owned ≈ expected
    real_state = [1.0,0.0]
    promoted = BlueTangle._variational_gate!(real_state,Op("Y",1))
    @test promoted == [0,im]
    @test real_state == [1.0,0.0]
    @test promoted !== real_state

    # Differentiate real and imaginary amplitudes, not only a norm-preserving
    # scalar. Finite differences use the old expanded operator as a reference.
    state = normalize(randn(rng, ComplexF64, 16))
    cases = [(Op("U3",_U3,3;control=1), [0.21,0.32,0.43]),
             (Op("FSIM",_FSIM,4,1;control=2), [0.27,0.38]),
             (adjoint(Op("U3",_U3,2)), [0.24,0.35,0.46])]
    for (op, parameters) in cases
        components(v) = [real.(v); imag.(v)]
        fast(p) = components(BlueTangle._variational_gate(state,op,p))
        reference(p) = components(op.expand(4,p...)*state)
        @test fast(parameters) ≈ reference(parameters)
        gradient = BlueTangle.ForwardDiff.jacobian(fast,parameters)
        finite_difference = hcat([begin
            step = zeros(length(parameters)); step[i] = 1e-6
            (reference(parameters+step)-reference(parameters-step))/2e-6
        end for i in eachindex(parameters)]...)
        @test gradient ≈ finite_difference atol=1e-9
    end

    # An encoding gate introduces Dual amplitudes; subsequent fixed gates must
    # preserve those derivatives without converting the state to ComplexF64.
    encoding = Op("U3",_U3,4)
    fixed = Op("CX",4,1;control=2)
    parameters = [0.21,0.32,0.43]
    encoded_loss(p) = begin
        encoded = BlueTangle._variational_gate(state,encoding,p)
        result = BlueTangle._variational_gate(encoded,fixed)
        abs2(result[1]) + 0.3abs2(result[12])
    end
    reference_loss(p) = begin
        result = fixed.expand(4)*encoding.expand(4,p...)*state
        abs2(result[1]) + 0.3abs2(result[12])
    end
    finite_difference = [begin
        step = zeros(3); step[i] = 1e-6
        (reference_loss(parameters+step)-reference_loss(parameters-step))/2e-6
    end for i in eachindex(parameters)]
    @test encoded_loss(parameters) ≈ reference_loss(parameters)
    @test BlueTangle.ForwardDiff.gradient(encoded_loss,parameters) ≈ finite_difference atol=1e-9
    owned_loss(p) = begin
        result = copy(state)
        result = BlueTangle._variational_gate!(result,encoding,p)
        result = BlueTangle._variational_gate!(result,fixed)
        abs2(result[1]) + 0.3abs2(result[12])
    end
    @test owned_loss(parameters) ≈ reference_loss(parameters)
    @test BlueTangle.ForwardDiff.gradient(owned_loss,parameters) ≈ finite_difference atol=1e-9

    @test_throws DimensionMismatch BlueTangle._variational_gate(ones(3),Op("H",1))
    @test_throws ArgumentError BlueTangle._variational_gate(ones(2),Op("H",2))
    @test_throws ArgumentError BlueTangle._variational_gate(ones(2),Op("H",0))
    @test_throws ArgumentError BlueTangle._variational_gate(ones(2),Op("H",1;control=0))
    @test_throws ArgumentError BlueTangle._variational_gate(ones(2),Op("MZ",1))
    @test_throws DimensionMismatch BlueTangle._variational_gate(ones(4),Op("bad",theta -> ones(3,3),1,2),[0.1])
end
