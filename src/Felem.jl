module Felem

include("GaussQuadrature.jl")
include("ShapeFunct.jl")

using .GaussQuadrature
using .ShapeFunct

export felem

function felem(q::Function, Le::Float64, nnel::Int, xstart::Float64=0.0; ngauss::Int=nnel)
    pem = ShapeFunctions(nnel)
    points, weights = gauss_quadrature(ngauss)
    J = Le / 2
    n_shape = length(pem)
    F_vect = zeros(n_shape)

    dof_scale = [isodd(i) ? 1.0 : J for i in 1:n_shape]

    for i in 1:n_shape
        val = 0.0
        for k in eachindex(points)
            x_actual = xstart + J * (1 + points[k])
            val += q(x_actual) * weights[k] * pem[i](points[k])
        end
        F_vect[i] = val * dof_scale[i] * J
    end
    
    return F_vect
end

end