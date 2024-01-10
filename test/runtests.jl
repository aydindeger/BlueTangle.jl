using BlueTangle
using Test

@testset "BlueTangle.jl" begin
    @test BlueTangle.int2bit(2,4)==[0,0,1,0]
end
