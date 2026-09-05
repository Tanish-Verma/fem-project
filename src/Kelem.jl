module Kelem

include("GaussQuadrature.jl")
include("ShapeFunct.jl")

using .GaussQuadrature
using .ShapeFunct

export kelem

function kelem(E::Function, I::Function, Le::Float64, nnel::Int, xstart::Float64=0.0; ngauss::Int=2*nnel-2)

    pem = ShapeFunctions(nnel)
    pemdashdash = (pem')'
    J = Le / 2
    n_shape = length(pem)
    K_mat = zeros(n_shape, n_shape)
    
    points, weights = gauss_quadrature(ngauss)
    
    dof_scale = [isodd(i) ? 1.0 : J for i in 1:n_shape]
    
    for i in 1:n_shape
        for j in 1:n_shape
            val = 0.0
            for k in eachindex(points)
                x_actual = xstart + J * (1 + points[k])
                
                N_ii = pemdashdash[i](points[k])
                N_jj = pemdashdash[j](points[k])
                
                val += weights[k] * E(x_actual) * I(x_actual) * N_ii * N_jj
            end
            K_mat[i, j] = val * dof_scale[i] * dof_scale[j]
        end
    end
    
    return K_mat * (1 / (J^3))
end

end