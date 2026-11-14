"""
    steane_prepare(basis::Symbol=:zero; blocks=1, backend=:statevector,
                   noise=false, prep_error=0.0, measurement_error=0.0,
                   fault=nothing, max_rounds=6, max_attempts=20)

Prepare independent logical |0⟩ or |+⟩ Steane blocks using physical product
preparation followed by verified Shor syndrome extraction and recovery.
Unlike `code.encode`, this path uses no general encoding circuit. Starting
from |0⟩ on every wire fixes logical Z=+1; the X checks project onto random
syndrome sectors, which Z recovery brings to the code space. For `:plus`,
physical Hadamards prepare |+⟩ on every wire and interchange X and Z roles.

`noise`, `prep_error`, `measurement_error`, and `fault` have the conventions
of `shor_qec_round`. Preparation and the optional physical H gates appear in
the same trace as extraction, with stages `:data_prep` and `:data_basis`.
The return value is the Shor round result with an additional `basis` field.
Choose `backend=:tensor` for a native MPS trajectory.

Under the single-fault model, an output may retain one correctable physical
error. This prepares stabilizer states only; it is not a magic-state factory.
"""
function steane_prepare(basis::Symbol=:zero; blocks::Int=1,
                        backend::Symbol=:statevector,
                        noise::Union{Bool,NoiseModel}=false,
                        prep_error::Real=0.0, measurement_error::Real=0.0,
                        fault=nothing, max_rounds::Int=6, max_attempts::Int=20)
    basis in (:zero, :plus) || throw(ArgumentError("basis must be :zero or :plus"))
    backend in (:statevector, :tensor) || throw(ArgumentError("backend must be :statevector or :tensor"))
    0 <= prep_error <= 1 || throw(ArgumentError("prep_error must lie between zero and one"))
    0 <= measurement_error <= 1 || throw(ArgumentError("measurement_error must lie between zero and one"))
    code = steane_code(blocks)
    state = zero_state(backend == :statevector ? code.n : N_MPS(code.n))
    trace = NamedTuple[]
    context = (stage=:data_prep, round=0, check=0, attempt=0)
    for q in 1:code.n
        state, _ = _qec_run_location(state, Op("I", q; noisy=false), context, trace;
                                     kind=:prep, prep_error=prep_error, fault=fault)
    end
    if basis == :plus
        context = (stage=:data_basis, round=0, check=0, attempt=0)
        for q in 1:code.n
            state, _ = _qec_run_location(state, Op("H", q), context, trace;
                                         noise=noise, fault=fault)
        end
    end
    result = _shor_qec_round(code, state, trace; noise=noise, prep_error=prep_error,
                             measurement_error=measurement_error, fault=fault,
                             max_rounds=max_rounds, max_attempts=max_attempts)
    return merge(result, (basis=basis,))
end
