using Test

include(joinpath(@__DIR__, "..", "src", "PostProcessing.jl"))

using .PostProcessing

@testset "Postprocessing of a simple beam" begin
    result = (
        u = zeros(4),
        mesh = (nodeLocs = [0.0, 1.0], nElem = 1),
        LM = reshape([1, 2, 3, 4], 4, 1),
        reactions = [10.0, 0.0, 0.0, 0.0],
    )
    bc = (
        nBC = 1,
        bcLoc = [0.0],
        bcDOF = [:w],
    )
    fm = (
        q = _ -> 2.0,
        pfLoc = Float64[],
        pfVal = Float64[],
    )

    results = postprocess(result, nothing, bc, fm; npoints = 5)

    @test length(results.xdef) == 50
    @test length(results.wdef) == 50
    @test first(results.xdef) == 0.0
    @test last(results.xdef) == 1.0
    @test results.wdef == zeros(50)
    @test results.x == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test results.V == [-10.0, -10.5, -11.0, -11.5, -12.0]
    @test results.M ≈ [0.0, -2.5625, -5.25, -8.0625, -11.0]
    @test results.rloc == [0.0]
    @test results.rval == [-10.0]
end