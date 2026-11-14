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
    coarse = apply(hamiltonian_exp(N,0.23,terms;dt=0.05),zero_state(N))
    fine = apply(hamiltonian_exp(N,0.23,terms;dt=0.005),zero_state(N))
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
