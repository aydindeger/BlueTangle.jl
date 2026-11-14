# OpenQASM import and export

`from_qasm` converts an OpenQASM 2 program into a vector of BlueTangle
operations. You can inspect or edit that vector, apply it to a state, and choose
a simulation backend. `to_qasm` writes the supported operations back to text or a
file.

## Import a string

```@example qasm
using BlueTangle

source = """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0], q[1];
"""
ops = from_qasm(source; is_file=false)
state = apply(zero_state(get_N(ops)), ops)
@assert correlation(state, "X,X", [1, 2]) ≈ 1
@assert correlation(state, "Z,Z", [1, 2]) ≈ 1
(qubits=get_N(ops), gates=length(ops))
```

OpenQASM indices start at 0, whereas BlueTangle indices start at 1. The importer
translates between these conventions. BlueTangle's statevector basis places
physical qubit 1 in the most significant bit.

## Read a file and round-trip a circuit

```@example qasm
restored = mktempdir() do directory
    path = joinpath(directory, "bell.qasm")
    to_qasm(ops, path)
    from_qasm(path; is_file=true)
end
@assert apply(zero_state(get_N(restored)), restored) ≈ state
println(to_qasm(restored))
```

For an existing file, use `from_qasm("circuits/bell.qasm"; is_file=true)`.
Without `is_file`, an existing path is read as a file and recognizable program
text is parsed directly. Explicit `is_file=false` is useful when a string should
always be interpreted as source.

## Registers and broadcasts

Multiple quantum registers are combined in declaration order. Declared qubits
with no gate operations are preserved: a four-qubit register remains four
qubits even if gates act only on `q[0]`. A trailing identity may be added to
represent that register size in the operation list.

Whole-register operations expand in register-index order. Equal-sized registers
can be paired, and a scalar operand can be broadcast against a register:

```@example qasm_registers
using BlueTangle
ops = from_qasm("""
OPENQASM 2.0;
qreg data[2];
qreg anc[2];
h data;
cx data, anc;
""")
state = apply(zero_state(get_N(ops)), ops)
@assert get_N(ops) == 4
@assert correlation(state, "Z,Z", [1, 3]) ≈ 1
@assert correlation(state, "Z,Z", [2, 4]) ≈ 1
(data_register=[1,2], anc_register=[3,4])
```

Here the `data` register maps to BlueTangle qubits 1–2 and `anc` maps to qubits
3–4. Operands can be full register names or individual indexed elements;
register slices are unsupported. Register sizes and indices are checked.

## Measurement records and reset

An imported `measure` becomes an explicit mid-circuit `MZ`.
`run(state, ops)` returns the final state and those outcomes in operation order:

```@example qasm_measurements
using BlueTangle
source = """
OPENQASM 2.0;
qreg q[2];
creg c[2];
x q[1];
measure q -> c;
reset q;
"""
ops = from_qasm(source)
result = run(zero_state(2), ops)
@assert result.record == [0, 1]
@assert result.state ≈ zero_state(2)

# The batch form returns one record per shot, starting from the zero state.
records = run(source, 3)
@assert records == [[0,1],[0,1],[0,1]]
records
```

The destination classical bits are validated during import. Their assignments
are not retained in the operation list; `run` reports measurement outcomes in
operation order. Classical `if` statements are unsupported.
To construct feed-forward circuits in Julia, use BlueTangle's conditional `Op`
or `ifOp`.

Use `sample(ops, shots)` for raw final outcomes or `measure(ops, shots)` for their
statistics. Both accept `backend=:statevector`, `:density_matrix`, or `:tensor`.
In the example above, `sample(ops, 3)` returns `[0, 0, 0]`: reset leaves the final
state in `|00⟩`, although the earlier record is `[0,1]`.

## Numeric parameters and diagnostics

Parameters accept finite real arithmetic, `pi`/`π`, implicit multiplication, and
supported numeric functions such as `sin`, `cos`, `sqrt`, `exp`, and `ln`:

```@example qasm_parameters
using BlueTangle
ops = from_qasm("qreg q[1]; ry(sin(pi/4) + 2^-2) q[0];")
@assert only(ops).mat ≈ _RY(sin(pi/4) + 0.25)
```

These expressions are parsed without executing Julia code. Assignments, macros,
module access, indexing, and arbitrary function calls are not accepted.

Errors identify the source, line, and statement:

```@example qasm_diagnostic
using BlueTangle
try
    from_qasm("OPENQASM 2.0;\nqreg q[2];\nx q[2];")
catch error
    println(sprint(showerror, error))
end
```

## Supported subset

The importer supports `OPENQASM 2.0`, `include "qelib1.inc"`, `qreg`/`creg`,
barriers, reset, and measurement, with these gates:

| Family | Names |
|:--|:--|
| One qubit | `id`, `x`, `y`, `z`, `h`, `s`, `sdg`, `t`, `tdg`, `sx`, `sxdg` |
| Parameterized one qubit | `rx`, `ry`, `rz`, `p`, `u1`, `u2`, `u3`, `U` |
| Controlled | `cx`/`cnot`, `cy`, `cz`, `ch`, `crx`, `cry`, `crz`, `cp`/`cu1`, `cu2`, `cu3` |
| Two qubits | `swap`, `iswap`, `fswap`, `syc`, `ecr`, `rxx`, `ryy`, `rzz`, `rxy`, `givens`, `fsim`, `swapa` |
| Three qubits | `ccx`, `ccy`, `ccz`, `cswap` |

Custom gate definitions, classical conditions, arbitrary include files, and
OpenQASM 3 constructs require a separate conversion before import. Unsupported
statements raise parsing errors.

Export preserves the represented gate matrices within floating-point tolerance.
For X/Y measurements, it writes the basis rotation, Z measurement, and inverse
rotation, so subsequent quantum operations see the same post-measurement state.
Some nonstandard gates are declared `opaque`; another tool must know their
semantics to execute them. In particular, BlueTangle uses
`RXX(theta) = exp(-im*theta*X⊗X)` and the analogous full-angle `RYY` convention;
check angle conventions when exchanging those gates with another simulator.

## Reference

```@docs
from_qasm
to_qasm
```
