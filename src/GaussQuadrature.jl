module GaussQuadrature

include("Polynomial.jl")
using .Polynomial
using LinearAlgebra

export gauss_quadrature, legendre_polynomial

function legendre_polynomial(n::Int)
    n >= 0 || throw(ArgumentError("n must be >= 0"))

    P0 = Poly([1.0])
    n == 0 && return P0

    xpoly = Poly([0.0, 1.0])
    P1 = xpoly
    n == 1 && return P1

    Pprev, Pcur = P0, P1
    for k in 1:(n-1)
        Pnext = ((2k+1)/(k+1)) * (xpoly * Pcur) - (k/(k+1)) * Pprev
        Pprev, Pcur = Pcur, Pnext
    end
    return Pcur
end

function gauss_quadrature(n::Int)
    n >= 1 || throw(ArgumentError("n must be >= 1"))

    Pn = legendre_polynomial(n)
    coeffs = Pn.elements
    lead = coeffs[end]
    a = zeros(n)
    for i in 1:n
        a[i] = coeffs[i] / lead
    end

    C = zeros(n, n)
    for i in 1:(n-1)
        C[i+1, i] = 1.0
    end
    C[:, n] = -a

    points = sort(real(eigvals(C)))

    weights = [2.0 / ((1 - x^2) * Pn'(x)^2) for x in points]

    return points, weights
end

end