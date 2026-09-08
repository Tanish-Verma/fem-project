using Test

include(joinpath(@__DIR__, "..", "src", "preprocessor.jl"))

using  .preprocessor

@testset "Mesh generation and LM matrix" begin
	beam = (nElem = 4,)
	bc = (bcLoc = [0.0, 4.0, 8.0],)
	fm = (pfLoc = Float64[], pmLoc = Float64[])
	rel = (relLoc = [4.0], relType = [:m])

	mesh = generateMesh(beam,bc, fm, rel;);
	LM = generateLM(mesh)

	println("Generated node locations: ", mesh.nodeLocs)
	println("Generated LM matrix:\n", LM)

	@test mesh.nodeLocs == [0.0, 2.0, 4.0, 6.0, 8.0]
	@test mesh.nElem == 4
	@test mesh.releaseElemIdx == [3]
	@test mesh.releaseNodeType == [:m]

	@test LM == [
		1 3 5 8;
		2 4 7 9;
		3 5 8 10;
		4 6 9 11;
	]
end

@testset "LM matrix for three-node elements" begin
	beam = (nElem = 4,)
	bc = (bcLoc = [0.0, 4.0, 8.0],)
	fm = (pfLoc = Float64[], pmLoc = Float64[])
	rel = (relLoc = [4.0], relType = [:m])

	mesh = generateMesh(beam, bc, fm, rel)
	LM = generateLM(mesh, 3)

	println("Three-node LM matrix:\n", LM)

	@test size(LM) == (6, mesh.nElem)
	@test LM == [
		1 5 9 14;
		2 6 11 15;
		3 7 12 16;
		4 8 13 17;
		5 9 14 18;
		6 10 15 19;
	]
end


