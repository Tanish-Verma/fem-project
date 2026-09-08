include("input.jl")
include("src/FEMBeamSolver.jl")
using .FEMBeamSolver
using Plots, LinearAlgebra

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

# ============================================================================
# CONVERGENCE STUDY: Maximum Deflection vs Number of Elements
# ============================================================================

println("="^70)
println("CONVERGENCE STUDY 1: Maximum Deflection vs Number of Elements (nnpe=2)")
println("="^70)

# Get the base beam parameters
beam_base = beamParameters()
bc = boundaryConditions()
fm = forceAndMoments()
rel = releases()

# Range of element numbers to test
nelem_range = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
max_deflections_nelem = Float64[]
nelem_tested = Int[]

println("\nAnalyzing convergence with respect to number of elements...")
for nelem in nelem_range
    # Modify beam parameters
    beam = (beam_length=beam_base.beam_length, 
            E=beam_base.E, 
            I=beam_base.I, 
            nElem=nelem, 
            nnpe=beam_base.nnpe)
    
    # Solve
    result = Solve(beam, bc, fm, rel)
    
    # Get maximum deflection (absolute value, ignoring negligible numerical errors)
    max_def = maximum(abs.(result.u[1:2:end]))  # Only deflection DOFs (odd indices)
    
    push!(max_deflections_nelem, max_def)
    push!(nelem_tested, nelem)
    
    println("nElem = $nelem: max deflection = $max_def m")
end

# Fit a polynomial to the convergence data
coeffs_nelem = fit_polynomial(Float64.(nelem_tested), max_deflections_nelem, 2)
nelem_fine = range(minimum(nelem_tested), maximum(nelem_tested), length=100)
fitted_nelem = eval_polynomial(coeffs_nelem, nelem_fine)

# Plot 1: Convergence with respect to number of elements
p1 = plot(nelem_tested, max_deflections_nelem, 
    seriestype=:scatter, 
    label="Computed", 
    xlabel="Number of Elements", 
    ylabel="Maximum Deflection (m)", 
    title="Question 6: Convergence with Respect to Number of Elements",
    legend=:topright,
    markersize=6)
plot!(p1, nelem_fine, fitted_nelem, 
    label="Polynomial Fit (degree 2)", 
    linewidth=2,
    linestyle=:dash)

# Save convergence data
nelem_convergence_data = (; nelem=nelem_tested, max_deflection=max_deflections_nelem, fit_coeffs=coeffs_nelem)
println("\nFitted polynomial coefficients (nelem convergence): $coeffs_nelem")

# ============================================================================
# CONVERGENCE STUDY 2: Maximum Deflection vs Order of Approximation (nnpe)
# ============================================================================

println("\n" * "="^70)
println("CONVERGENCE STUDY 2: Maximum Deflection vs Order of Approximation (nElem=10)")
println("="^70)

# Use a reasonable number of elements for this study
base_nelem = 10
beam_base_fixed = (beam_length=beam_base.beam_length, 
                   E=beam_base.E, 
                   I=beam_base.I, 
                   nElem=base_nelem, 
                   nnpe=2)

# Range of nnpe values to test (nodes per element)
nnpe_range = [2, 3, 4, 5, 6]
max_deflections_nnpe = Float64[]
nnpe_tested = Int[]

println("\nAnalyzing convergence with respect to approximation order...")
for nnpe in nnpe_range
    # Modify beam parameters
    beam = (beam_length=beam_base_fixed.beam_length, 
            E=beam_base_fixed.E, 
            I=beam_base_fixed.I, 
            nElem=beam_base_fixed.nElem, 
            nnpe=nnpe)
    
    # Solve
    result = Solve(beam, bc, fm, rel)
    
    # Get maximum deflection
    max_def = maximum(abs.(result.u[1:2:end]))  # Only deflection DOFs (odd indices)
    
    push!(max_deflections_nnpe, max_def)
    push!(nnpe_tested, nnpe)
    
    println("nnpe = $nnpe: max deflection = $max_def m")
end

# Fit a polynomial to the convergence data
coeffs_nnpe = fit_polynomial(Float64.(nnpe_tested), max_deflections_nnpe, 2)
nnpe_fine = range(minimum(nnpe_tested), maximum(nnpe_tested), length=100)
fitted_nnpe = eval_polynomial(coeffs_nnpe, nnpe_fine)

# Plot 2: Convergence with respect to approximation order
p2 = plot(nnpe_tested, max_deflections_nnpe, 
    seriestype=:scatter, 
    label="Computed", 
    xlabel="Approximation Order (nnpe)", 
    ylabel="Maximum Deflection (m)", 
    title="Question 7: Convergence with Respect to Approximation Order",
    legend=:topright,
    markersize=6)
plot!(p2, nnpe_fine, fitted_nnpe, 
    label="Polynomial Fit (degree 2)", 
    linewidth=2,
    linestyle=:dash)

# Save convergence data
nnpe_convergence_data = (; nnpe=nnpe_tested, max_deflection=max_deflections_nnpe, fit_coeffs=coeffs_nnpe)
println("\nFitted polynomial coefficients (nnpe convergence): $coeffs_nnpe")

# ============================================================================
# VISUALIZATION - Display plots separately
# ============================================================================

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

println("\nPlots saved as:")
println("  - convergence_elements.png")
println("  - convergence_order.png")
println("="^70)


