@testset "Numeric gate expressions and QASM" begin
    for (text, value) in (("2pi/3", 2pi/3), (".25π", pi/4), ("sin(pi/2)+cos(0)", 2.0),
                          ("sqrt(2)^2", 2.0), ("2(1+3)/4", 2.0), ("-exp(0)+ln(e)", 0.0), ("2^-2", 0.25))
        @test gates("RX($text)") ≈ _RX(value)
    end
    @test gates("U3(sin(pi/2),sqrt(2),atan(1,2))") ≈ _U3(1.0, sqrt(2), atan(1,2))
    for text in ("begin; 1; end", "(x=1)", "Base.pi", "getfield(Base,:pi)",
                 "[1][1]", "error(\"executed\")", "@eval 1", "true", "Inf", "1/0")
        @test_throws ArgumentError gates("RX($text)")
    end
    @test_throws ArgumentError gates("RX()")
    @test_throws ArgumentError gates("RX(1,2)")
    @test_throws ArgumentError gates("unknown")
    @test all(op -> op.q == 1 && op.qubit == 1, random_clifford(1,8))
    @test BlueTangle.random_gate(1).qubit == 1

    names1 = ["I", "X", "Y", "Z", "H", "S", "SD", "T", "TD", "SX", "SXD",
              "RX(sin(pi/3))", "RY(0.2)", "RZ(2pi)", "P(0.3)", "U1(0.3)",
              "U2(0.2,0.4)", "U3(0.1,0.2,0.3)"]
    names2 = ["CX", "CY", "CZ", "SWAP", "ISWAP", "FSWAP", "SYC", "ECR",
              "CP(0.2)", "RXX(0.3)", "RYY(0.3)", "RZZ(0.3)", "RXY(0.3)",
              "GIVENS(0.3)", "FSIM(0.2,0.3)", "SWAPA(0.2)"]
    roundtrip_ops = QuantumOps[Op(name,1) for name in names1]
    append!(roundtrip_ops, [Op(name,1,3) for name in names2])
    append!(roundtrip_ops, [Op(name,1;control=3) for name in ["X", "Y", "Z", "H", "S", "SD", "T", "TD", "RX(0.3)", "RY(0.3)", "RZ(0.3)", "U1(0.3)", "U2(0.2,0.3)", "U3(0.1,0.2,0.3)"]])
    append!(roundtrip_ops, [Op(name,1,2,3) for name in ("CCX", "CCY", "CCZ", "CSWAP")])
    initial = normalize(ComplexF64[1, 2im, 3, -1im, 2, 1+im, -2, 4im])
    for op in roundtrip_ops
        imported = from_qasm(to_qasm([op]; nqubits=3))
        @test get_N(imported) == 3
        @test apply(initial, imported) ≈ apply(initial, op)
    end

    # Non-Z measurements must preserve the post-measurement state too.
    for name in ("MX", "M(X)", "MY", "M(Y)", "MZ", "M(Z)")
        ops = [Op(name,1), Op("H",1)]
        Random.seed!(314)
        direct, values = apply(initial, ops; track_measurements=true)
        Random.seed!(314)
        imported, qasm_values = apply(initial, from_qasm(to_qasm(ops; nqubits=3)); track_measurements=true)
        @test values == qasm_values
        @test direct ≈ imported
    end

    multireg = from_qasm("OPENQASM 2.0; qreg data[2]; qreg anc[1]; x anc[0];")
    @test get_N(multireg) == 3
    @test apply(zero_state(3), multireg) ≈ apply(zero_state(3), Op("X",3))
    @test get_N(from_qasm("OPENQASM 2.0; qreg q[4]; x q[0];")) == 4
    @test get_N(from_qasm("OPENQASM 2.0; qreg q[2];")) == 2
    @test_throws ArgumentError from_qasm("qreg q[2]; x q[2];")
    @test_throws ArgumentError from_qasm("qreg q[2]; x other[0];")
    @test_throws ArgumentError from_qasm("gate fake a { x a; } qreg q[1]; fake q[0];")
    @test_throws ArgumentError from_qasm("rx(error(\"executed\")) q[0];")
    @test_throws ArgumentError to_qasm([Op("RZZ(0.2)",1,2;control=3)])
    @test_throws ArgumentError to_qasm([Op("X",gate.Z,1)])
    @test_throws ArgumentError to_qasm([Op("X",2)]; nqubits=1)

    broadcast = from_qasm("qreg a[2]; qreg b[2]; h a; cx a,b;";is_file=false)
    expected_ops = [Op("H",1),Op("H",2),Op("CX",1,3),Op("CX",2,4)]
    @test apply(zero_state(4), broadcast) ≈ apply(zero_state(4), expected_ops)
    measurement = from_qasm("qreg q[2]; creg c[2]; x q; measure q -> c; reset q;")
    measured_state, record = apply(zero_state(2), measurement;track_measurements=true)
    @test record == [1,1]
    @test measured_state ≈ zero_state(2)
    @test only(from_qasm("U(pi,0,pi) q[0];")).mat ≈ gate.X
    @test_throws ArgumentError from_qasm("qreg q[2]; creg c[1]; measure q -> c;")
    @test_throws ArgumentError from_qasm("qreg q[2]; creg c[2]; measure q[0] -> c[2];")
    @test_throws ArgumentError from_qasm("qreg q[2]; h q";is_file=false)
    bad_source = "OPENQASM 2.0;\n// keep source line numbers\nqreg q[2];\nx q[2];\n"
    error = try from_qasm(bad_source;is_file=false); catch e; e; end
    @test error isa ArgumentError
    @test occursin("<string>:4:",sprint(showerror,error))
    @test occursin("x q[2]",sprint(showerror,error))
    mktempdir() do directory
        path = joinpath(directory,"circuit.qasm")
        to_qasm([Op("H",1),Op("CX",1,2)],path)
        @test apply(zero_state(2), from_qasm(path;is_file=true)) ≈ apply(zero_state(2), [Op("H",1),Op("CX",1,2)])
        write(path,bad_source)
        error = try from_qasm(path); catch e; e; end
        @test occursin("$path:4:",sprint(showerror,error))
    end
end
