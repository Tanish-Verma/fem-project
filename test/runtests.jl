using Test

include(joinpath(@__DIR__, "..", "src", "Kelem.jl"))
include(joinpath(@__DIR__, "..", "src", "Felem.jl"))
include(joinpath(@__DIR__, "testpreprocessor.jl"))

using .Kelem: kelem
using .Felem: felem

function analytical_stiffness(EI, Le)
    return (EI / Le^3) .* [
         12       6*Le   -12       6*Le
          6*Le    4*Le^2  -6*Le    2*Le^2
        -12      -6*Le    12      -6*Le
          6*Le    2*Le^2  -6*Le    4*Le^2
    ]
end

function analytical_uniform_load(q, Le)
    return (q * Le / 2) .* [1, Le/6, 1, -Le/6]
end

@testset "Two-node Euler-Bernoulli element" begin
    E0 = 210.0e9
    I0 = 100.0e-6
    q0 = -90.0e3

    for Le in (2.0, 4.0)
        K = kelem(_ -> E0, _ -> I0, Le, 2)
        F = felem(_ -> q0, Le, 2)

        @test K ≈ analytical_stiffness(E0 * I0, Le) rtol=1e-12
        @test K ≈ transpose(K) rtol=1e-14
        @test F ≈ analytical_uniform_load(q0, Le) rtol=1e-12
    end
end

@testset "Element position and higher order" begin
    E(x) = 2.0 + x
    I(x) = 3.0

    K_at_zero = kelem(E, I, 2.0, 2, 0.0; ngauss=3)
    K_shifted = kelem(E, I, 2.0, 2, 1.0; ngauss=3)
    K_higher = kelem(_ -> 1.0, _ -> 1.0, 3.0, 3)
    F_higher = felem(_ -> 2.0, 3.0, 3)

    @test K_at_zero != K_shifted
    @test size(K_higher) == (6, 6)
    @test K_higher ≈ transpose(K_higher) atol=1e-13
    @test length(F_higher) == 6
    @test sum(F_higher[1:2:end]) ≈ 6.0 atol=1e-12
end
