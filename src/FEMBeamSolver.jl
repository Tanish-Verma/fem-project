module FEMBeamSolver

# Include all supporting modules in dependency order
include("Polynomial.jl")
include("GaussQuadrature.jl")
include("ShapeFunct.jl")
include("Kelem.jl")
include("Felem.jl")
include("Global_KF.jl")
include("preprocessor.jl")
include("Solver.jl")
include("PostProcessing.jl")

# Use all submodules to bring their exports into this namespace
using .Polynomial
using .GaussQuadrature
using .ShapeFunct
using .Kelem
using .Felem
using .Global_KF
using .preprocessor
using .Solver
using .PostProcessing

# Export main solving functions
export Solve, postprocess, plot_results, generateMesh, generateLM, global_kf, kelem, felem, gauss_quadrature, ShapeFunctions

end # module FEMBeamSolver
