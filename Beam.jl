include("input.jl")
include("src/FEMBeamSolver.jl")
using .FEMBeamSolver

beam = beamParameters()
bc = boundaryConditions()
fm = forceAndMoments()
rel = releases()

result = Solve(beam,bc,fm,rel)

println("Displacement vector:")
println(result.u)
println("\nReaction vector:")
println(result.reactions)

results = postprocess(result,beam,bc,fm)

plots = plot_results(results)
display(plots.combined_plot)
println("Press Enter to close the plot...")
readline()
