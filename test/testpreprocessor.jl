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

@testset "LM matrix for three-node elements" begin
	beam = beamParameters()
	bc = boundaryConditions()
	fm = forceAndMoments()
	rel = releases()

	mesh = generateMesh(beam, bc, fm, rel)
	LM = generateLM(mesh, 3)

	println("Three-node LM matrix:\n", LM)

	@test size(LM) == (6, mesh.nElem)
	@test LM == [
		1 5 10 14;
		2 7 11 15;
		3 8 12 16;
		4 9 13 17;
		5 10 14 18;
		6 11 15 19;
	]
end


