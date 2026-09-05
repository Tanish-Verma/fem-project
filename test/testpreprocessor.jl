using Test

include(joinpath(@__DIR__, "..", "src", "preprocessor.jl"))
include(joinpath(@__DIR__,"..","input.jl"));

using  .preprocessor

@testset "Mesh generation and LM matrix" begin
	beam = beamParameters()
	bc = boundaryConditions()
	fm = forceAndMoments()
	rel = releases()

	mesh = generateMesh(beam,bc, fm, rel;);
	LM = generateLM(mesh)

	println("Generated node locations: ", mesh.nodeLocs)
	println("Generated LM matrix:\n", LM)

	@test mesh.nodeLocs == [0.0, 2.0, 4.0, 6.0, 8.0]
	@test mesh.nElem == 4
	@test mesh.releaseElemIdx == [2]
	@test mesh.releaseNodeType == [:m]

	@test LM == [
		1 3 6 8;
		2 5 7 9;
		3 6 8 10;
		4 7 9 11;
	]
end


