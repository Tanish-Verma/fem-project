include("input.jl")
include("src/Solver.jl")
include("src/PostProcessing.jl")
using .Solver
using .PostProcessing

beam = beamParameters()
bc = boundaryConditions()
fm = forceAndMoments()
rel = releases()

result = Solve(beam,bc,fm,rel)
results = postprocess(result,beam,bc,fm)

println("\nDisplacement vector:")
println(result.u)
println("\nReaction vector:")
println(result.reactions)

plots = plot_results(results,result)
