using BlueTangle
using Test
using LinearAlgebra

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

    @test begin #test sa.sparse of a circuit

        N=10
        depth=50
        ops=random_ops(N,depth)
        c=compile(ops)

        #let's test this on a random state
        s=random_state(N)
        s1=deepcopy(s)
        s2=deepcopy(s)
        for o=ops
            s1=apply(s1,o)
        end

        U=sa.sparse(c)
        s2=U*s2

        isapprox(s1,s2,atol=1e-12)
    end

    @test begin # test decomposition
        
        U = gate.H
        α, β, γ, δ = zyz_decomposition(U)
        U2 = exp(im*α) * _RZ(β) * _RY(γ) * _RZ(δ)
        t1=isapprox(la.norm(U - U2),0,atol=1e-10)

        C = kron(gates("RZ(.1pi)"), gates("RY(.3)"))
        A, B = kronecker_decomposition(C)
        t2=isapprox(la.norm(kron(A, B) - C),0,atol=1e-10)
        t1 && t2

    end

    @test begin # density matrix

        N=4
        depth=20
        ops=random_ops(N,depth);
        nm=NoiseModel("depolarizing",.1)

        #exact
        sym="T"

        state=zero_state(N)
        for o=ops
            state=apply(state,o)
        end
        mag_e=sum(expect(state,sym))

        rho=zero_state(N)*zero_state(N)'
        for o=ops
            rho=apply(rho,o)
        end
        mag_r=sum(expect(state,sym))

        b1=isapprox(mag_e,mag_r,atol=1e-10);
        b2=isapprox(state*state',rho,atol=1e-10);


        mag_list=[]
        exp_no=1000
        for experiment=1:exp_no
            state_n=zero_state(N)
            for o=ops
                state_n=apply(state_n,o;noise=nm)
            end
            push!(mag_list,sum(expect(state_n,sym)))
        end

        rho_n=zero_state(N)*zero_state(N)'
        for o=ops
            rho_n=apply(rho_n,o;noise=nm)
        end

        b3=isapprox(sum(expect(rho_n,sym)),mean(mag_list),atol=10/exp_no)

        b1 && b2 && b3
    end

end
