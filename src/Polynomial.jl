module Polynomial

export Poly, derivative

struct Poly
    elements::Vector{Float64}
end

import Base: +, -, *, /, adjoint, zero, transpose

Base.zero(::Type{Poly}) = Poly([0.0])
Base.zero(p::Poly) = Poly([0.0])
function -(p::Poly)
    return Poly(-p.elements)
end

function +(p1::Poly, p2::Poly)
    deg = length(p1.elements) > length(p2.elements) ? length(p1.elements) : length(p2.elements)
    p = Poly(zeros(deg))
    for i in 1:deg
        p.elements[i] = i<=length(p1.elements) ? p1.elements[i] : 0
        p.elements[i] += i<=length(p2.elements) ? p2.elements[i] : 0
    end
    return p
end

-(p1::Poly, p2::Poly) = p1 + (-p2)

*(c::Number, p::Poly) = Poly(c.*p.elements)

*(p::Poly, c::Number) = c*p

/(p::Poly, c::Number) = p*(1/c)

function *(p1::Poly, p2::Poly)
    deg = length(p1.elements) + length(p2.elements) - 2
    p = Poly(zeros(deg + 1))
    for i in 1:length(p1.elements)
        for j in 1:length(p2.elements)
            p.elements[i + j - 1] += p1.elements[i] * p2.elements[j]
        end
    end
    return p
end

function (p::Poly)(x::Number)
    result = 0.0
    for i in length(p.elements):-1:1
        result = p.elements[i] + result * x
    end
    return result
end

function derivative(p::Poly)
    if length(p.elements) == 1
        return Poly([0.0])
    end
    deg = length(p.elements) - 1
    p_dash = Poly(zeros(deg))
    for i in 2:length(p.elements)
        p_dash.elements[i-1] = (i-1) * p.elements[i]
    end
    return p_dash
end

Base.adjoint(p::Poly) = derivative(p)

Base.transpose(p::Poly) = p

end