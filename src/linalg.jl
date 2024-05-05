function la.ishermitian(op::Op)
    if op.mat isa AbstractMatrix && la.ishermitian(op.mat)
        return true
    else
        return false
    end
end
function Base.adjoint(op::Op)
    matf = op.mat
    name = op.name * "†"

    if matf isa Function
        adjmatf = x -> adjoint(matf(x))
    elseif matf isa la.Adjoint
        name = op.name[1:end-1]
        adjmatf = matf.parent
    else
        adjmatf = adjoint(matf)
    end
    return Op(name, adjmatf, op.qubit, op.target_qubit; type=op.type, noisy=op.noisy, control=op.control)
end

function Base.adjoint(circ::Circuit)
    all_layers = circ.layers
    adjoint_all_layers = []

    while !isempty(all_layers)
        layer = pop!(all_layers)
        adjointlayer = map(layer) do op
            return la.ishermitian(op) ? op : adjoint(op)
        end
        push!(adjoint_all_layers, adjointlayer)
    end

    adjointcirc = deepcopy(circ)

    copy!(adjointcirc.layers, adjoint_all_layers)

    return adjointcirc
end
