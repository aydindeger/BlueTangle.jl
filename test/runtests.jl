using BlueTangle
using Test

@testset "BlueTangle.jl" begin
    @test BlueTangle.int2bin(2,4)==[0,0,1,0]

    @test begin #this checks nonlocal hilbert construction
        N=6
        qubit=2
        control=5
        a=Op("CNOT",qubit,control).expand(N)
        b=Op("X",control;control=qubit).expand(N)
        a==b
    end

    @test begin #nonlocal hilbert + noise channel
        N=6
        mixed_state=sa.sparse(la.normalize(rand(Complex{Float64}, 2^N)));
        nm=QuantumChannel(2,"depolarizing",.5)
        kraus_vec=nm.kraus;
        
        qubit=2
        second_qubit=3
        
        pA=partial_trace(mixed_state,[qubit,second_qubit])
        probs=[real(la.tr(kraus * pA * kraus')) for kraus in kraus_vec]
        probs2=[sum(abs2.(hilbert(N,kraus,qubit,second_qubit)*mixed_state)) for kraus in kraus_vec]
        
        op=QC("my quantum channel",kraus_vec,qubit,second_qubit)
        probs3=op.prob(mixed_state)
        probs ≈ probs2 ≈ probs3
    end

end
