module PostProcessing

include("ShapeFunct.jl")
using .ShapeFunct
using Plots

gr()
# pythonplot()

export postprocess, plot_results

function postprocess(result, beam, bc, fm; npoints=200)
    u = result.u
    mesh = result.mesh
    LM = result.LM
    xdef = Float64[]
    wdef = Float64[]
    xtheta = Float64[]
    theta = Float64[]

    for e in 1:mesh.nElem
        x1 = mesh.nodeLocs[e]
        x2 = mesh.nodeLocs[e+1]
        Le = x2-x1
        J = Le/2
        lm = LM[:,e]
        ue = u[lm]
        nnpe = div(length(ue),2)
        N = ShapeFunctions(nnpe)
        Ndash = [p' for p in N]
        points = range(-1,1,length=50)

        for ξ in points
            x = x1 + J*(1+ξ)
            w = 0.0
            rot = 0.0
            for i in eachindex(N)
                if isodd(i)
                    w += N[i](ξ)*ue[i]
                    rot += Ndash[i](ξ)*ue[i]/J
                else
                    w += N[i](ξ)*J*ue[i]
                    rot += Ndash[i](ξ)*ue[i]
                end
            end
            push!(xdef,x)
            push!(wdef,w)
            push!(xtheta,x)
            push!(theta,rot)
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

    # --- Correct for the missing constant from clamped-end moment reactions ---
    for i in 1:bc.nBC
        if bc.bcDOF[i] != :θ
            continue
        end
        loc = bc.bcLoc[i]
        if isapprox(loc, xmin)
            node = findfirst(x -> isapprox(x,loc), mesh.nodeLocs)
            dof = LM[2,node]
            M .+= result.reactions[dof]   # note: NOT negated, opposite sign to :w reactions
        end
    end

    pmLoc = hasproperty(fm, :pmLoc) ? fm.pmLoc : Float64[]
    pmVal = hasproperty(fm, :pmVal) ? fm.pmVal : Float64[]
    for (loc,moment) in zip(pmLoc,pmVal)
        for i in eachindex(x)
            if x[i] >= loc
                M[i] += moment
            end
        end
    end

    return (; xdef,wdef,xtheta,theta,x,V,M,rloc,rval)
end


function plot_results(results)

    common = (linewidth=2, tickfontsize=9, guidefontsize=10,
              titlefontsize=12, legend=:best)

    x_closed = [results.x[1]; results.x; results.x[end]]
    V_closed = [0.0; results.V; 0.0]
    M_closed = [0.0; results.M; 0]

    p1 = plot(results.xdef, results.wdef;
        xlabel="Position x (m)", ylabel="Deflection (m)",
        title="Deflection Profile", label="Deflection",
        yformatter=:scientific, common...)

    ptheta = plot(results.xtheta, results.theta;
        xlabel="Position x (m)", ylabel="Rotation θ (rad)",
        title="Rotation Profile", label="Rotation θ",
        yformatter=:scientific, common...)

    p2 = plot(x_closed, V_closed;
        xlabel="Position x (m)", ylabel="Shear Force (kN)",
        title="Shear Force Diagram", label="Shear Force",
        fill=(0, 0.5), fillcolor=:pink,
        common...)

    p3 = plot(x_closed, M_closed;
        xlabel="Position x (m)", ylabel="Bending Moment (kN m)",
        title="Bending Moment Diagram", label="Bending Moment",
        fill=(0, 0.5), fillcolor=:green,
        common...)

    combined = plot(p1, ptheta, p2, p3,
        layout=(2,2),
        size=(1400,1000),
        left_margin=8Plots.mm,
        right_margin=5Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=8Plots.mm)

    savefig(combined,"combined_results.pdf")

    return (; deflection_plot=p1,
        rotation_plot=ptheta,
        shear_plot=p2,
        moment_plot=p3,
        combined_plot=combined)
end

end
