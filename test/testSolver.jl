using Test

include(joinpath(@__DIR__, "..", "src", "Solver.jl"))

using .Solver

@testset "One-element cantilever under uniform load" begin
	beam = (
		beam_length = 2.0,
		E = _ -> 1.0,
		I = _ -> 1.0,
		nElem = 1,
		nnpe = 2,
	)
	bc = (
		nBC = 2,
		bcLoc = [0.0, 0.0],
		bcDOF = [:w, :θ],
		bcVal = [0.0, 0.0],
	)
	fm = (
		q = _ -> 1.0,
		nPointForce = 1,
		nPointMoment = 0,
		pfLoc = [2.0],
		pfVal = [0.0],
		pmLoc = Float64[],
		pmVal = Float64[],
	)
	rel = (nReleases = 0, relLoc = Float64[], relType = Symbol[])

	sol = Solve(beam, bc, fm, rel)

	@test sol.mesh.nodeLocs == [0.0, 2.0]
	@test sol.LM == reshape([1, 2, 3, 4], 4, 1)
	@test sol.u[1:2] ≈ [0.0, 0.0] atol = 1e-10
	@test sol.u[3] ≈ 2.0 atol = 1e-10
	@test sol.u[4] ≈ 4.0 / 3.0 atol = 1e-10
end
