@testset "Dependency-free circuit and measurement display" begin
    function rendered(args...; kwargs...)
        io = IOBuffer()
        @test isnothing(plotq(args...; io=io, kwargs...))
        return String(take!(io))
    end
    text = rendered([Op("H", 1), Op("CX", 1, 2), Op("RY(0.123456)", 2)])
    @test occursin("q1:", first(split(text, '\n')))
    @test occursin("●", text) && occursin("⊕", text)
    @test occursin("RY(0.123456)", text)
    @test occursin("RZ(-0.08)", rendered(Op("RZ(-0.08000000000000002)",1)))
    @test occursin("custom(0.123456789)", rendered(Op("custom(0.123456789)",gate.H,1)))
    @test occursin("empty circuit", rendered(QuantumOps[]))
    @test occursin("empty circuit", rendered([]))
    @test occursin("αβγδ", rendered(Op("αβγδ", gate.H, 1)))
    @test occursin("custom block", rendered(OpF("custom block", identity)))
    @test occursin("0: I; 1: X@q2", rendered(ifOp("MZ", 1, Op("I", 1), Op("X", 2))))
    @test occursin("q2:", first(split(rendered([Op("CX",1,2)]; reverse_bits=true), '\n')))
    @test occursin("data:", rendered([Op("H",1)]; labels=["data"]))
    @test_throws ArgumentError plotq([Op("H",1)]; labels=["a","b"])
    @test_throws ArgumentError plotq([Op("H",1)]; backend=:unknown)

    # Displaying a function-valued gate must never evaluate it or consume RNG.
    calls = Ref(0)
    gate_function(a, b) = (calls[] += 1; _FSIM(a, b))
    functional = Op("custom fSim", gate_function, 1, 2)
    rendered(functional)
    @test calls[] == 0
    Random.seed!(123)
    expected = rand()
    Random.seed!(123)
    rendered(functional)
    @test rand() == expected
    @test occursin("000", rendered(measure(zero_state(3))))
    @test occursin("Other", rendered(measure(plus_state(4)); max_outcomes=3))
    @test occursin("Uniform", rendered([measure(plus_state(2)), measure(zero_state(2))], ["Uniform", "Zero"]))
    reference = measure(zero_state(1))
    signed = Measurement([0,1], [1.1,-0.1], reference.expect, reference.mag_moments,
                         reference.measurement_basis, reference.number_of_experiment,
                         "Signed estimate", 1, reference.density_matrix)
    @test occursin("Measurement estimates (quasiprobabilities)", rendered(signed))

    if isnothing(Base.get_extension(BlueTangle, :BlueTanglePyPlotExt))
        @test_throws ArgumentError plotq([Op("H",1)]; backend=:pyplot)
    end
end

if !isnothing(Base.get_extension(BlueTangle, :BlueTanglePyPlotExt))
    @testset "Optional graphical display" begin
        plots = Base.get_extension(BlueTangle, :BlueTanglePyPlotExt).PyPlot
        mktempdir() do directory
            evaluations = Ref(0)
            not_evaluated(a, b) = (evaluations[] += 1; _FSIM(a, b))
            cases = [QuantumOps[],
                     QuantumOps[OpF("block", identity)],
                     QuantumOps[Op("H",1), Op("CX",1,3), Op("SWAP",2,3), ifOp("MZ",2,"I","X")],
                     QuantumOps[Op("custom", not_evaluated, 1, 2)],
                     QuantumOps[OpQC("amplitude_damping",0.2,1), Op("RY(0.37)",3; control=2)]]
            for (i, ops) in enumerate(cases)
                figure = plotq(ops; backend=:pyplot, fold=2, show=false)
                file = joinpath(directory, "circuit_$i.png")
                @test savefigure(figure, file) == file
                @test filesize(file) > 1000
            end
            @test evaluations[] == 0
            figure = plotq([Op("H",1),Op("H",2)]; backend=:pyplot, show=false)
            qubit_labels = [text for text in figure.axes[1].texts if text.get_text() in ("q1","q2")]
            @test qubit_labels[1].get_position()[2] > qubit_labels[2].get_position()[2]
            figure = plotq([measure(plus_state(2)),measure(zero_state(2))], ["Uniform","Zero"];
                           backend=:pyplot, show=false)
            @test length(figure.axes[1].patches) == 8
            @test [tick.get_text() for tick in figure.axes[1].get_xticklabels()] == ["00","01","10","11"]
            figure = plotq(measure(plus_state(6)); backend=:pyplot, max_outcomes=4, show=false)
            @test figure.axes[1].get_xticklabels()[end].get_text() == "Other"
            @test sum(bar.get_height() for bar in figure.axes[1].patches) ≈ 1
            reference = measure(zero_state(1))
            signed = Measurement([0,1], [1.1,-0.1], reference.expect, reference.mag_moments,
                                 reference.measurement_basis, reference.number_of_experiment,
                                 "Signed estimate", 1, reference.density_matrix)
            figure = plotq(signed; backend=:pyplot, show=false)
            @test figure.axes[1].get_ylim()[1] < -0.1
            @test [bar.get_height() for bar in figure.axes[1].patches] ≈ [1.1,-0.1]
            @test figure.axes[1].get_ylabel() == "Quasiprobability"
            figure = plotq([measure(plus_state(3)),measure(zero_state(3))]; backend=:pyplot,show=false)
            @test isempty(figure.axes[1].texts)
            # Loading PyPlot must not change the default backend.
            io = IOBuffer()
            @test isnothing(plotq([Op("H",1)]; io=io))
            @test occursin("[H]", String(take!(io)))
        end
        plots.close("all")
    end
end
