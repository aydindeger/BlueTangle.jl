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
        random_state=sa.sparse(la.normalize(rand(Complex{Float64}, 2^N)));
        nm=QuantumChannel(2,"depolarizing",.5)
        kraus_vec=nm.kraus;
        
        qubit=2
        second_qubit=3
        
        pA=partial_trace(random_state,[qubit,second_qubit])
        probs=[real(la.tr(kraus * pA * kraus')) for kraus in kraus_vec]
        probs2=[sum(abs2.(hilbert(N,kraus,qubit,second_qubit)*random_state)) for kraus in kraus_vec]
        
        op=OpQC("my quantum channel",kraus_vec,qubit,second_qubit)
        probs3=op.prob(random_state)
        probs ≈ probs2 ≈ probs3
    end

    @test begin #test stinespring_dilation

    kraus_ops=OpQC("amplitude_damping",.2,2,3).kraus
    mat = stinespring_dilation(kraus_ops);
    
    N=2
    state=random_state(N)
    rho=state*state'
    zero=zero_state(N)   
    rho_zero=zero*zero'
    a=sum(kraus * rho * kraus' for kraus=kraus_ops)
    b=mat * kron(rho_zero,rho) * mat'
    c=partial_trace(sa.sparse(b),[4,4],[1])
    a==c

    end

end
