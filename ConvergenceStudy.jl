include("input.jl")
include("src/Solver.jl")
include("src/ConvergenceMetrics.jl")
using .Solver: Solve
using .ConvergenceMetrics
using Plots

gr()
default(linewidth=2, framestyle=:box, grid=true, size=(800, 500), margin=5Plots.mm)

beam_base = beamParameters()
bc, fm, rel = boundaryConditions(), forceAndMoments(), releases()

special_locations = sort(unique(vcat(bc.bcLoc, fm.pfLoc, fm.pmLoc, rel.relLoc, 0.0, beam_base.beam_length)))
base_nelem = length(special_locations) - 1

function run_study(beams)
    return [
        begin
            res = Solve(beam, bc, fm, rel)
            (; nelem=res.mesh.nElem, degree=2*beam.nnpe-1, max_def=sampled_max_deflection(res))
        end
        for beam in beams
    ]
end

# Q6: h-refinement (fixed degree, doubling number of elements)
h_beams = [merge(beam_base, (; nElem=base_nelem * (2^i))) for i in 0:5]
h_data = run_study(h_beams)

p1 = plot([d.nelem for d in h_data], [d.max_def for d in h_data]; 
    marker=:circle, xlabel="Total number of elements", ylabel="Maximum |deflection| (m)",
    title="Maximum Deflection vs. Total Number of Elements", yscale=:log10, label="Max Deflection")

# Q7: p-refinement (fixed mesh of base_nelem elements, increasing degree)
p_beams = [merge(beam_base, (; nElem=base_nelem, nnpe=p)) for p in 2:10]
p_data = run_study(p_beams)

p2 = plot([d.degree for d in p_data], [d.max_def for d in p_data]; 
    marker=:circle, xlabel="Polynomial degree p = 2nnpe-1", ylabel="Maximum |deflection| (m)",
    title="Maximum Deflection vs. Higher-Order Approximation", yscale=:log10, label="Max Deflection")

savefig(p1, "question6_convergence.png")
savefig(p2, "question7_convergence.png")

if stdin isa Base.TTY
    display(p1)
    println("Showing plot 1. Press Enter to view plot 2...")
    readline()
    display(p2)
    println("Showing plot 2. Press Enter to finish...")
    readline()
elseif isinteractive()
    display(p1)
    display(p2)
end