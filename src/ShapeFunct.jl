module ShapeFunct
include("Polynomial.jl")
using LinearAlgebra
using .Polynomial
export ShapeFunctions

function ShapeFunctions(num)
    C = zeros(2*num, 2*num)
    p = [Poly([j==i ? 1.0 : 0.0 for j in 1:2*num]) for i in 1:2*num]
    x = collect(range(-1.0, 1.0, length=num))
    for j in 1:2*num
        C[1:2:end, j] = p[j].(x)
        C[2:2:end, j] = p[j]'.(x)
    end
    return inv(C') * p
end

# function ShapeFunctions(num)
#     C = zeros(2*num, 2*num)
#     p = [Poly([1.0]) for _ in 1:2*num]
#     p_dash = [derivative(p[i]) for i in 1:2*num]
#     x = collect(range(-1.0, 1.0, length=num))
#     for i in 1:num
#         C[2*i-1, :] = [x[i]^(j-1) for j in 1:2*num]
#         C[2*i, :] = [j==1 ? 0 : (j-1)*x[i]^(j-2) for j in 1:2*num]
#     end
#     return inv(C')
# end
end