@testset "Variational setup and gradients" begin
    given = [0.31]
    opt = AnsatzOptions(N=1, ops=["RY"], loss=gate.Z, pars_initial=given,
                        model="descent", number_of_iterations=3, history=false)
    @test opt.pars_initial == given
    @test opt.pars_initial !== given
    @test variational_apply(given,opt) ≈ _RY(given[1]) * zero_state(1)
    energy, parameters = VQA(opt)
    @test energy ≈ cos(only(parameters))
    @test energy < cos(only(given))
    @test given == [0.31]
    @test opt.pars_initial == given
    @test_throws DimensionMismatch variational_apply([0.1,0.2],opt)
    @test_throws DimensionMismatch AnsatzOptions(N=2,ops=["RY"],loss=state->0.0,pars_initial=[0.1])
    @test_throws ArgumentError AnsatzOptions(N=1,ops=["RY"],loss=gate.Z,model="typo")

    prepared = apply(zero_state(1), [Op("X",1)])
    opt_init = AnsatzOptions(N=1,ops=["RY"],loss=gate.Z,init=prepared,pars_initial=[0.0])
    @test variational_apply([0.0],opt_init) ≈ one_state(1)
    swap_options = AnsatzOptions(N=2,ops=["SWAPA"],loss=hilbert(2,gate.Z,1),pars_initial=[0.3])
    @test swap_options.dim == 1
    @test swap_options.ops[1].q == 2
    @test_throws ArgumentError hamiltonian(2,[1.0,"X",2.0])
    @test_throws ArgumentError hamiltonian(2,[1.0,"X"],"typo")
    @test_throws ArgumentError hamiltonian((0,2),[1.0,"X"])
    adjoint_options = AnsatzOptions(N=1,ops=[adjoint(Op("U3",_U3,1))],loss=gate.Z,pars_initial=[0.2,0.3,0.4])
    @test adjoint_options.dim == 3
    @test variational_apply([0.2,0.3,0.4],adjoint_options) ≈ _U3(0.2,0.3,0.4)' * zero_state(1)
    @test BlueTangle.loss_and_grad_paramshift([0.2,0.3,0.4],adjoint_options)[2] ≈ [-sin(0.2),0,0] atol=1e-10

    # Independent finite differences check every built-in parameter-shift rule,
    # including gates with a different angle scale and controlled rotations.
    for (name, fn) in BlueTangle._parameterized_gates
        count = BlueTangle._find_argument_number(fn)
        p = collect(range(0.21, step=0.17, length=count))
        two_qubit = size(fn(p...),1) == 4
        for controlled in (false,true)
            n = two_qubit ? 3 : 2
            op = two_qubit ? Op(name,fn,1,2;control=controlled ? 3 : -2) :
                             Op(name,fn,1;control=controlled ? 2 : -2)
            init = normalize(ComplexF64[complex(i,(-1)^i) for i in 1:2^n])
            observable = hilbert(n,gate.X,1) + 0.3hilbert(n,gate.Z,n)
            options = AnsatzOptions(N=n,ops=[op],loss=observable,init=init,pars_initial=p,paramshift=true)
            _, gradient = BlueTangle.loss_and_grad_paramshift(p,options)
            finite_difference = [begin
                step = zeros(length(p)); step[i] = 1e-6
                (options.loss(variational_apply(p+step,options))-options.loss(variational_apply(p-step,options)))/2e-6
            end for i in eachindex(p)]
            @test gradient ≈ finite_difference atol=1e-8
        end
    end

    for model in ("adam", "lbfgs")
        options = AnsatzOptions(N=1,ops=["RY"],loss=gate.Z,pars_initial=[0.3],
                                model=model,number_of_iterations=model == "lbfgs" ? 12 : 3,history=false,paramshift=true)
        energy, p = VQA(options)
        @test energy ≈ cos(only(p))
    end
end
