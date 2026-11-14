@testset "Operation dependencies and layers" begin
    conditional = ifOp("MZ", 1, Op("I", 3), Op("X", 3))
    ops = [Op("X", 1), conditional, Op("H", 3)]
    direct = apply(ops, zero_state(3))
    @test expect(direct, Op("X", 3)) ≈ -1
    @test apply(zero_state(3), ops) ≈ direct
    @test get_N([conditional]) == 3
    @test get_stats([conditional]).N == 3
    @test get_stats(ops).depth == length(get_layers(ops)) == 3

    # A branch's control and target qubits both constrain subsequent operations.
    controlled_branch = ifOp("MZ", 1, Op("X", 2; control=4), Op("CX", 3, 5))
    @test get_N([controlled_branch]) == 5
    @test length(get_layers([Op("H", 4), controlled_branch, Op("Z", 5)])) == 3

    controlled_ops = [Op("X", 3), Op("X", 1; control=3), Op("H", 1)]
    @test get_stats(controlled_ops).depth == length(get_layers(controlled_ops)) == 3
    @test apply(zero_state(3), controlled_ops) ≈ apply(controlled_ops, zero_state(3))

    # Independent operations can still share a layer; arbitrary functions cannot.
    parallel_ops = [Op("X", 1), Op("X", 2), Op("CX", 1, 2), Op("Z", 3)]
    @test get_stats(parallel_ops).depth == length(get_layers(parallel_ops)) == 2
    barrier = OpF("identity", identity)
    barrier_ops = [Op("X", 1), barrier, Op("H", 2)]
    @test get_stats(barrier_ops).depth == length(get_layers(barrier_ops)) == 3
    @test length(get_layers([barrier])) == 1
    @test apply(zero_state(2), barrier_ops) ≈ apply(barrier_ops, zero_state(2))

    @test isempty(get_layers(Op[]))
    @test get_stats(Op[]).depth == 0
    @test isempty(get_optimized_layers([Op("X", 1), Op("X", 1)]))
end
