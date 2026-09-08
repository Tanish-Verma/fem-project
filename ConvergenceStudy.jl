include("input.jl")
include("src/FEMBeamSolver.jl")
using .FEMBeamSolver
using Plots, LinearAlgebra

gr()  # ensure GR backend (crisp, good for saving)
default(
    fontfamily = "Computer Modern",
    titlefontsize = 14,
    guidefontsize = 12,
    tickfontsize = 10,
    legendfontsize = 10,
    linewidth = 2,
    framestyle = :box,
    grid = true,
    gridalpha = 0.3,
    size = (800, 550),
    margin = 5Plots.mm,
)

# Function to fit a polynomial of degree d to data
function fit_polynomial(x::Vector, y::Vector, d::Int)
    A = hcat([x.^i for i in 0:d]...)
    coeffs = A \ y
    return coeffs
end

# Function to evaluate fitted polynomial
function eval_polynomial(coeffs::Vector, x)
    return sum(coeffs[i] .* x.^(i-1) for i in 1:length(coeffs))
end


beam_base = beamParameters()
bc = boundaryConditions()
fm = forceAndMoments()
rel = releases()

nelem_range = [5, 10, 15, 20, 30, 50, 100]
max_deflections_nelem = Float64[]
nelem_tested = Int[]

for nelem in nelem_range
    # Modify beam parameters
    beam = (beam_length=beam_base.beam_length, 
            E=beam_base.E, 
            I=beam_base.I, 
            nElem=nelem, 
            nnpe=beam_base.nnpe)
    
    # Solve
    result = Solve(beam, bc, fm, rel)
    
    max_def = maximum(abs.(result.u[1:2:end]))
    
    push!(max_deflections_nelem, max_def)
    push!(nelem_tested, nelem)
    
    println("nElem = $nelem: max deflection = $max_def m")
end

# Fit a polynomial to the convergence data
coeffs_nelem = fit_polynomial(Float64.(nelem_tested), max_deflections_nelem, 2)
nelem_fine = range(minimum(nelem_tested), maximum(nelem_tested), length=100)
fitted_nelem = eval_polynomial(coeffs_nelem, nelem_fine)


p1 = plot(nelem_tested, max_deflections_nelem, 
    seriestype = :scatter, 
    label = "Computed", 
    xlabel = "Number of Elements", 
    ylabel = "Maximum Deflection (m)", 
    title = "Question 6: Convergence with Respect to Number of Elements",
    legend = :topleft,
    markersize = 7,
    markerstrokewidth = 0.5,
    markercolor = :dodgerblue,
)
plot!(p1, nelem_fine, fitted_nelem, 
    label = "Polynomial Fit (degree 2)", 
    linewidth = 2.5,
    linestyle = :dash,
    linecolor = :crimson,
)

nelem_convergence_data = (; nelem=nelem_tested, max_deflection=max_deflections_nelem, fit_coeffs=coeffs_nelem)

base_nelem = 5
beam_base_fixed = (beam_length=beam_base.beam_length, 
                   E=beam_base.E, 
                   I=beam_base.I, 
                   nElem=base_nelem, 
                   nnpe=2)

# Range of nnpe values to test (nodes per element)
nnpe_range = [2, 3, 4, 5, 6, 7, 8, 9, 10]
max_deflections_nnpe = Float64[]
nnpe_tested = Int[]
for nnpe in nnpe_range
    beam = (beam_length=beam_base_fixed.beam_length, 
            E=beam_base_fixed.E, 
            I=beam_base_fixed.I, 
            nElem=beam_base_fixed.nElem, 
            nnpe=nnpe)
    
    result = Solve(beam, bc, fm, rel)
    
    max_def = maximum(abs.(result.u[1:2:end]))
    
    push!(max_deflections_nnpe, max_def)
    push!(nnpe_tested, nnpe)
    
    println("nnpe = $nnpe: max deflection = $max_def m")
end

coeffs_nnpe = fit_polynomial(Float64.(nnpe_tested), max_deflections_nnpe, 2)
nnpe_fine = range(minimum(nnpe_tested), maximum(nnpe_tested), length=100)
fitted_nnpe = eval_polynomial(coeffs_nnpe, nnpe_fine)

p2 = plot(nnpe_tested, max_deflections_nnpe, 
    seriestype = :scatter, 
    label = "Computed", 
    xlabel = "Approximation Order (nnpe)", 
    ylabel = "Maximum Deflection (m)", 
    title = "Question 7: Convergence with Respect to Approximation Order",
    legend = :topleft,
    markersize = 7,
    markerstrokewidth = 0.5,
    markercolor = :dodgerblue,
)
plot!(p2, nnpe_fine, fitted_nnpe, 
    label = "Polynomial Fit (degree 2)", 
    linewidth = 2.5,
    linestyle = :dash,
    linecolor = :crimson,
)

nnpe_convergence_data = (; nnpe=nnpe_tested, max_deflection=max_deflections_nnpe, fit_coeffs=coeffs_nnpe)
println("\nFitted polynomial coefficients (nnpe convergence): $coeffs_nnpe")

# Save high-resolution figures to disk
savefig(p1, "question6_convergence.png")
savefig(p2, "question7_convergence.png")

# Display plot 1
display(p1)
println("\nPress Enter to close the plots...")
readline()

# Display plot 2
display(p2)
println("\nPress Enter to close the plots...")
readline()

println("\n" * "="^70)
println("SUMMARY")
println("="^70)
println("\nQuestion 6 - Element Convergence:")
println("  Range: $nelem_range elements")
println("  Max deflection range: $(minimum(max_deflections_nelem)) to $(maximum(max_deflections_nelem)) m")
println("  Convergence behavior: Polynomial fit degree 2")

println("\nQuestion 7 - Order Convergence:")
println("  Range: $nnpe_range nodes per element")
println("  Max deflection range: $(minimum(max_deflections_nnpe)) to $(maximum(max_deflections_nnpe)) m")
println("  Convergence behavior: Polynomial fit degree 2")

println("="^70)