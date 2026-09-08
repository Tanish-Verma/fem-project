module PostProcessing

include("ShapeFunct.jl")
using .ShapeFunct
using Plots
gr()

export postprocess, plot_results

function postprocess(result, beam, bc, fm; npoints=200)
    u = result.u
    mesh = result.mesh
    LM = result.LM
    xdef = Float64[]
    wdef = Float64[]

    for e in 1:mesh.nElem
        x1 = mesh.nodeLocs[e]
        x2 = mesh.nodeLocs[e+1]
        Le = x2-x1
        J = Le/2
        lm = LM[:,e]
        ue = u[lm]
        nnpe = div(length(ue),2)
        N = ShapeFunctions(nnpe)
        points = range(-1,1,length=50)

        for ξ in points
            x = x1 + J*(1+ξ)
            w = 0.0
            for i in eachindex(N)
                if isodd(i)
                    w += N[i](ξ)*ue[i]
                else
                    w += N[i](ξ)*J*ue[i]
                end
            end
            push!(xdef,x)
            push!(wdef,w)
        end
    end

    rloc = Float64[]
    rval = Float64[]

    for i in 1:bc.nBC
        if bc.bcDOF[i] != :w
            continue
        end
        loc = bc.bcLoc[i]
        node = findfirst(x -> isapprox(x,loc),mesh.nodeLocs)
        if node < length(mesh.nodeLocs)
            dof = LM[1,node]
        else
            dof = LM[end-1,end]
        end
        push!(rloc,loc)
        push!(rval,-result.reactions[dof])
    end

    xmin = mesh.nodeLocs[1]
    xmax = mesh.nodeLocs[end]
    x = collect(range(xmin,xmax,length=npoints))
    V = zeros(npoints)

    # Compute cumulative integral of distributed load using trapezoidal rule
    q_integral = zeros(npoints)
    for i in 2:npoints
        dx = x[i] - x[i-1]
        q_integral[i] = q_integral[i-1] + (fm.q(x[i-1]) + fm.q(x[i])) * dx / 2
    end

    for i in eachindex(x)
        xi = x[i]
        for j in eachindex(rloc)
            if rloc[j] <= xi
                V[i] += rval[j]
            end
        end

        V[i] -= q_integral[i]

        for j in eachindex(fm.pfLoc)
            if fm.pfLoc[j] <= xi
                V[i] -= fm.pfVal[j]
            end
        end
    end

    M = zeros(npoints)

    for i in 2:npoints
        dx = x[i]-x[i-1]
        M[i] = M[i-1] + (V[i-1]+V[i])*dx/2
    end
    return (; xdef,wdef,x,V,M,rloc,rval)
end


function plot_results(results)

    p1 = plot(results.xdef,results.wdef,
        xlabel="Position x (m)",
        ylabel="Deflection (m)",
        title="Deflection Profile",
        label="Deflection")

    p2 = plot(results.x,results.V,
        xlabel="Position x (m)",
        ylabel="Shear Force (kN)",
        title="Shear Force Diagram",
        label="Shear Force")

    p3 = plot(results.x,results.M,
        xlabel="Position x (m)",
        ylabel="Bending Moment (kN m)",
        title="Bending Moment Diagram",
        label="Bending Moment")

    combined = plot(p1,p2,p3,
        layout=(3,1),
        size=(800,1000))

    return (; deflection_plot=p1,
        shear_plot=p2,
        moment_plot=p3,
        combined_plot=combined)
end

end
