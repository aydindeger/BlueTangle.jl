@testset "Backend integration and MPS ordering" begin
    rng = MersenneTwister(92)
    sites = N_MPS(4)
    vector = normalize(randn(rng, ComplexF64, 16))
    psi = to_MPS(vector, sites; cutoff=0.0)
    @test to_state(psi) ≈ vector
    @test inner(psi, im*vector) ≈ im
    @test inner(im*vector, psi) ≈ -im
    @test fidelity(psi,psi) ≈ 1
    @test correlation(psi,"X,Y",[1,4]) ≈ correlation(vector,"X,Y",[1,4])
    @test expect(psi,Op("Z",2)) ≈ expect(vector,Op("Z",2))
    before = to_state(psi)
    @test length(sample_bit(psi, 5)) == 5
    @test to_state(psi) ≈ before
    @test entanglement_entropy(zero_state(N_MPS(1))) == 0
    @test sample_exact(psi)[2] ≈ abs2.(vector)
    @test measure(psi).number_of_experiment == 1000

    large = one_state(N_MPS(72))
    measurement = measure(large, 3)
    @test measurement.bitstr == [big(2)^72-1]
    @test measurement.expect == -ones(72)
    @test size(measurement.density_matrix) == (0,0)
    @test_throws ArgumentError sample_exact(large)
    @test int2bin(bin2int(ones(Int,72)),72) == ones(Int,72)

    for backend in (:statevector, :tensor, :density_matrix)
        circuit = compile([Op("I",1)],Options(readout_noise=Noise1("bit_flip",1.0)))
        @test measure(circuit, 20; backend=backend).bitstr == [1]
        @test run([Op("X",1),Op("MZ",1),Op("RES",1),Op("MZ",1)],2;backend=backend) == [[1,0],[1,0]]
        reset = compile([Op("H",1),Op("CX",1,2),Op("RES",1)], Options(density_matrix=true))
        result = measure(reset, 1000;backend=backend)
        @test result.expect[1] ≈ 1
        @test abs(result.expect[2]) < 0.15
        @test tr(result.density_matrix) ≈ 1
    end
    # A Kraus channel must be sampled afresh on each trajectory, even when no
    # NoiseModel option or measurement gate is present in the circuit.
    noisy = compile([Op("I",1),OpQC("bit_flip",0.5,1)])
    @test abs(measure(noisy,1000).expect[1]) < 0.15
    @test abs(measure(noisy,1000;backend=:density_matrix).expect[1]) < 0.15

    mixed = ComplexF64[0.25 0; 0 0.75]
    @test measure(mixed).sample ≈ [0.25,0.75]
    @test measure(mixed).expect ≈ [-0.5]
    @test (Op("X",1) * mixed) ≈ gate.X * mixed * gate.X

    options = Options(readout_noise=Noise1("bit_flip",0.1),measurement_mitigate=true)
    mitigated = measure(compile([Op("I",1)],options),1)
    @test length(mitigated.sample) == 2
    @test sum(mitigated.sample) ≈ 1
    @test minimum(mitigated.sample) < 0
    @test abs(only(mitigated.expect)) ≈ 1.25
    # A finite-shot histogram need not lie in the image of physical input
    # probabilities. Inversion must retain signed weights to undo the response.
    _, signed = BlueTangle._mitigate_probabilities([0],[1.0],1,options.readout_noise)
    response = [0.9 0.1; 0.1 0.9]
    @test signed ≈ [1.125,-0.125]
    @test response * signed ≈ [1.0,0.0] atol=1e-14
    _, signed_pair = BlueTangle._mitigate_probabilities([0],[1.0],2,options.readout_noise)
    @test kron(response,response) * signed_pair ≈ [1.0,0.0,0.0,0.0] atol=1e-14
    bits, corrected = BlueTangle._mitigate_probabilities([0,1,2,3],[0.81,0.09,0.09,0.01],2,options.readout_noise)
    @test bits == 0:3
    @test corrected ≈ [1,0,0,0] atol=1e-14

    for op in (Op("custom",gate.CX,1,3), Op("custom",gate.CY,3,1), Op("SWAP",1,3;control=2))
        @test sa.sparse(compile([op])) ≈ op.expand(3)
    end
    conditional = ifOp("MZ",1,Op("custom",gate.CX,1,3),Op("I",1))
    ops = [Op("H",3),conditional]
    @test to_state(compile(ops)) ≈ apply(ops,zero_state(3))
    @test isempty([name for name in names(BlueTangle) if !isdefined(BlueTangle,name)])
end

@testset "Common operation constructors" begin
    conditional = Op("MZ",1;if0="I",if1="X")
    @test conditional isa ifOp
    @test apply(one_state(1),conditional) ≈ zero_state(1)
    @test Op(Noise1("bit_flip",1.0),1) isa OpQC
    @test apply(zero_state(1),Op(Noise1("bit_flip",1.0),1)) ≈ one_state(1)
    block = Op("prepare",[Op("X",3)])
    @test block isa OpF
    @test get_N([block]) == 3
    @test to_state(compile([block])) ≈ apply(zero_state(3),Op("X",3))
    @test apply(zero_state(1),Op("full gate",gate.X)) ≈ one_state(1)
    @test apply(zero_state(1),Op("identity transform",identity)) ≈ zero_state(1)
    @test_throws ArgumentError Op("X",1;if1="X")
    @test !isdefined(BlueTangle,:measure_ZNE)
end

@testset "Callback sampling and controlled routing" begin
    calls = Ref(0)
    callback = Op("alternating preparation", state -> begin
        calls[] += 1
        isodd(calls[]) ? apply(state, Op("X",1)) : state
    end)
    circuit = compile([Op("I",1), callback])
    result = measure(circuit,10;backend=:density_matrix)
    @test calls[] == 10
    @test result.expect == [0.0]
    op = Op("X",1;control=4)
    routed = compile([op];layout=Layout(ones(Int,1,4)))
    @test sa.sparse(routed) ≈ op.expand(4)
    @test all(g -> g isa Op && (g.q == 1 ? abs(g.control-g.qubit)==1 : abs(g.qubit-g.target_qubit)==1), Iterators.flatten(routed.layers))
end
