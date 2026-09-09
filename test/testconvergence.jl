using Test
include("../src/ConvergenceMetrics.jl")
include("../src/Solver.jl")
include("../input.jl")
using .ConvergenceMetrics
using .Solver: Solve

@testset "Convergence metrics" begin
    # w(x) = x - x^2/2 on [0,2]: zero nodal deflections, interior max 1/2.
    result = (u=[0.0, 1.0, 0.0, -1.0], LM=reshape(1:4, 4, 1),
        mesh=(nElem=1, nodeLocs=[0.0, 2.0]))
    @test sampled_max_deflection(result; samples=10000) ≈ 0.5 atol=1e-4
    @test_throws ArgumentError sampled_max_deflection(result; samples=1)
end

@testset "h and p refinement" begin
    beam = beamParameters()
    bc, fm, rel = boundaryConditions(), forceAndMoments(), releases()
    coarse = Solve(merge(beam, (; nElem=3, nnpe=2)), bc, fm, rel)
    fine = Solve(merge(beam, (; nElem=6, nnpe=2)), bc, fm, rel)
    higher = Solve(merge(beam, (; nElem=3, nnpe=3)), bc, fm, rel)
    @test fine.mesh.nodeLocs[1:2:end] ≈ coarse.mesh.nodeLocs
    @test maximum(diff(fine.mesh.nodeLocs)) ≈ maximum(diff(coarse.mesh.nodeLocs))/2
    @test higher.mesh.nodeLocs == coarse.mesh.nodeLocs
    @test all(r -> isfinite(sampled_max_deflection(r)), (coarse, fine, higher))
end