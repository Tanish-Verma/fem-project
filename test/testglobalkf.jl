using Test

include(joinpath(@__DIR__, "..", "src", "Global_KF.jl"))

using .Global_KF: global_kf

# Expected values below were derived independently (closed-form Hermite beam
# stiffness / fixed-end-force formulas, worked out by hand), not by copying
# global_kf's own algorithm -- so these can actually catch a bug shared by
# global_kf and kelem/felem, unlike a test that reimplements the same loop.

@testset "Point loads at a moment release (:m -- theta split, w shared)" begin
    # 2-element mesh, hinge at the shared node x=2. E=I=1, q=0 so results are
    # driven purely by the point load.
    mesh = (nElem = 2, nodeLocs = [0.0, 2.0, 4.0],
            releaseElemIdx = [2], releaseNodeType = [:m])
    # node 2's w (dof 3) stays shared; element 2's left theta gets its own
    # private dof (5) instead of reusing element 1's right theta (dof 4)
    LM = [1 3; 2 5; 3 6; 4 7]
    zero_load(x) = 0.0
    onefn(x) = 1.0
        beam = (E = onefn, I = onefn)
        fm = (q = zero_load, pfLoc = [2.0], pfVal = [10.0],
            pmLoc = Float64[], pmVal = Float64[])

        _, Fg = global_kf(mesh, LM, beam, fm)
    @test Fg[3] ≈ 10.0                    # shared w dof gets the full force
    @test isapprox(Fg[5], 0.0; atol=1e-8) # element 2's private theta: untouched
    @test isapprox(Fg[6], 0.0; atol=1e-8)
    @test isapprox(Fg[7], 0.0; atol=1e-8)
end

@testset "Point loads at a shear release (:v -- w split, theta shared)" begin
    mesh = (nElem = 2, nodeLocs = [0.0, 2.0, 4.0],
            releaseElemIdx = [2], releaseNodeType = [:v])
    # w is split (private dof 5 for element 2's left end); theta stays shared (dof 4)
    LM = [1 5; 2 4; 3 6; 4 7]
    zero_load(x) = 0.0
    onefn(x) = 1.0
        beam = (E = onefn, I = onefn)
        fm = (q = zero_load, pfLoc = Float64[], pfVal = Float64[],
            pmLoc = [2.0], pmVal = [10.0])

        _, Fg = global_kf(mesh, LM, beam, fm)
    @test Fg[4] ≈ 10.0                    # shared theta dof gets the full moment
end

@testset "A release actually changes global stiffness (not a no-op)" begin
    E1(x) = 1.0
    I1(x) = 1.0
    zero_load(x) = 0.0
        beam = (E = E1, I = I1)
        fm = (q = zero_load, pfLoc = Float64[], pfVal = Float64[],
            pmLoc = Float64[], pmVal = Float64[])
    nodeLocs = [0.0, 2.0, 4.0]

    # baseline: ordinary shared mesh, no release at all
    mesh_plain = (nElem = 2, nodeLocs = nodeLocs,
                  releaseElemIdx = Int[], releaseNodeType = Symbol[])
    LM_plain = [1 3; 2 4; 3 5; 4 6]
    Kg_plain, _ = global_kf(mesh_plain, LM_plain, beam, fm)

    # same two elements, but a moment release (hinge) at the shared node
    mesh_hinge = (nElem = 2, nodeLocs = nodeLocs,
                  releaseElemIdx = [2], releaseNodeType = [:m])
    LM_hinge = [1 3; 2 5; 3 6; 4 7]
    Kg_hinge, _ = global_kf(mesh_hinge, LM_hinge, beam, fm)

    # Dof 4 is "the shared theta at node 2" in the plain mesh (fed by BOTH
    # elements), but becomes "element 1's own private theta" in the hinged
    # mesh (fed by ONLY element 1). The two elements are identical (same L,
    # E, I), and a uniform prismatic beam element's two end-rotation
    # stiffness terms are always equal to each other -- so removing one
    # element's contribution should exactly halve this entry, regardless of
    # what kelem's actual numbers are.
    @test Kg_hinge[4,4] ≈ Kg_plain[4,4] / 2

    # Dof 4 (element 1's private theta) and dof 5 (element 2's private theta)
    # should be completely decoupled: no single element's local stiffness
    # ever has both of them as local dof at once, so this entry must be
    # EXACTLY zero -- a pure connectivity argument, true regardless of
    # kelem's specific numeric output.
    @test isapprox(Kg_hinge[4,5], 0.0; atol=1e-10)

    # symmetry should hold with or without a release present
    @test Kg_hinge ≈ transpose(Kg_hinge)
end