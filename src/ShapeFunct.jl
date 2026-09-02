module ShapeFunct
include("Polynomial.jl")
using LinearAlgebra
using .Polynomial
export ShapeFunctions, ShapeFunctionsPoly

function ShapeFunctionsPoly(num)
    C = zeros(2*num, 2*num)
    p = [Polynom(2*num-1, [j==i ? 1.0 : 0.0 for j in 1:2*num]) for i in 1:2*num]
    p_dash = [derivative(p[i]) for i in 1:2*num]
    x = collect(range(-1.0, 1.0, length=num))
    for i in 1:num
        C[2*i-1, :] = [eval(p[j], x[i]) for j in 1:2*num]
        C[2*i, :] = [eval(p_dash[j], x[i]) for j in 1:2*num]
    end
    C_ = inv(C')
    C_p = [Polynom(2*num-1, C_[i, :]) for i in 1:2*num]
    return C_p
end

function ShapeFunctions(num)
    C = zeros(2*num, 2*num)
    p = ones(2*num)
    p_dash = [j-1 for j in 1:2*num]
    x = collect(range(-1.0, 1.0, length=num))
    for i in 1:num
        C[2*i-1, :] = [x[i]^(j-1) for j in 1:2*num]
        C[2*i, :] = [j==1 ? 0 : (j-1)*x[i]^(j-2) for j in 1:2*num]
    end
    return inv(C')
end
end