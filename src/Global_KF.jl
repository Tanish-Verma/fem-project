module Global_KF

include("Kelem.jl")
include("Felem.jl")
include("ShapeFunct.jl")

using .Kelem
using .Felem
using .ShapeFunct

export global_kf

function global_kf(mesh, LM, beam, fm; nnpe::Int = div(size(LM, 1), 2))

    nElem = mesh.nElem
    nodeLocs = mesh.nodeLocs

    totdofpElem = size(LM, 1)
    if totdofpElem != 2 * nnpe
        error("Mismatch: LM matrix has $totdofpElem rows, but nnpe = $nnpe expects $(2*nnpe) rows.")
    end

    ndof = maximum(LM)

    Kg = zeros(ndof, ndof)
    Fg = zeros(ndof)

    # 1. ELEMENT ASSEMBLY (Stiffness & Continuous Load)
    for e in 1:nElem
        xstart = nodeLocs[e]
        xend   = nodeLocs[e+1]
        Le     = xend - xstart

        Ke = kelem(beam.E, beam.I, Le, nnpe, xstart)
        Fe = felem(fm.q, Le, nnpe, xstart)

        lm = LM[:, e]

        Kg[lm, lm] .+= Ke
        Fg[lm]     .+= Fe
    end

    # Helper function: maps point coordinate xp to element index
    function find_elem(xp)
        if isapprox(xp, nodeLocs[end]; atol=1e-9)
            return nElem
        end
        idx = findfirst(i -> nodeLocs[i] <= xp < nodeLocs[i+1], 1:nElem)
        if isnothing(idx) && isapprox(xp, nodeLocs[1]; atol=1e-9)
            return 1
        end
        return idx
    end

    # 2. POINT FORCES
    for (xp, P) in zip(fm.pfLoc, fm.pfVal)
        e = find_elem(xp)
        isnothing(e) && error("Point force at x = $xp is outside the mesh domain.")

        xstart, xend = nodeLocs[e], nodeLocs[e+1]
        Le = xend - xstart
        J  = Le / 2.0
        ξ  = 2.0 * (xp - xstart) / Le - 1.0

        pem = ShapeFunctions(nnpe)
        lm  = LM[:, e]

        # Rotation DOFs (even indices) are scaled by the Jacobian J -- same
        # convention kelem/felem use, needed here for exactly the same reason
        # it's needed in the point_moments loop below: pem[i] for a theta-type
        # shape function is calibrated to unit d/dxi, not unit d/dx.
        dof_scale = [isodd(i) ? 1.0 : J for i in 1:length(lm)]
        Fe_point = [dof_scale[i] * pem[i](ξ) * P for i in 1:length(lm)]
        Fg[lm] .+= Fe_point
    end

    # 3. POINT MOMENTS
    for (xp, M) in zip(fm.pmLoc, fm.pmVal)
        e = find_elem(xp)
        isnothing(e) && error("Point moment at x = $xp is outside the mesh domain.")

        xstart, xend = nodeLocs[e], nodeLocs[e+1]
        Le = xend - xstart
        J  = Le / 2.0
        ξ  = 2.0 * (xp - xstart) / Le - 1.0

        pem = ShapeFunctions(nnpe)

        # Polynomial derivative operator '
        pem_dash = [p' for p in pem]

        lm  = LM[:, e]

        # Rotation DOFs (even indices) are scaled by Jacobian J
        dof_scale = [isodd(i) ? 1.0 : J for i in 1:length(lm)]
        Fe_moment = [(dof_scale[i] / J) * pem_dash[i](ξ) * M for i in 1:length(lm)]
        Fg[lm] .+= Fe_moment
    end

    return Kg, Fg
end

end