module ConvergenceMetrics
include("ShapeFunct.jl")
using .ShapeFunct
export sampled_max_deflection

function sampled_max_deflection(result; samples=200)
    samples >= 2 || throw(ArgumentError("at least two samples required"))
    N = ShapeFunctions(size(result.LM, 1) ÷ 2)
    basis = [p(xi) for xi in range(-1.0, 1.0, length=samples), p in N]
    
    max_val = 0.0
    for e in 1:result.mesh.nElem
        J = (result.mesh.nodeLocs[e+1] - result.mesh.nodeLocs[e]) / 2
        ue = copy(result.u[result.LM[:, e]])
        ue[2:2:end] .*= J
        max_val = max(max_val, maximum(abs, basis * ue))
    end
    return max_val
end

end