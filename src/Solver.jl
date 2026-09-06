module Solver

include("Global_KF.jl")
include("preprocessor.jl")

using .preprocessor
using .Global_KF
using LinearAlgebra

export Solve

function Solve(beam, bc, fm, rel)
    # 1. Generate Mesh and Location Matrix
    mesh = generateMesh(beam, bc, fm, rel)
    LM = generateLM(mesh,beam.nnpe)
    
    # 2. Assemble Global K and F (assuming global_kf returns standard matrices)
    K, F = global_kf(mesh, LM, beam, fm, nnpe=beam.nnpe)
    K_orig = copy(K);
    F_orig = copy(F);
    # 3. Apply Boundary Conditions
    for i in 1:bc.nBC
        loc = bc.bcLoc[i]
        dof_type = bc.bcDOF[i]
        val = bc.bcVal[i]
        
        # A. Find which element this node belongs to
        node_idx = findfirst(==(loc), mesh.nodeLocs)
        
        # B. Grab the exact global DOF number from the LM
        global_dof = 0
        if node_idx < length(mesh.nodeLocs)
            # It's the LEFT node of element `node_idx`
            if dof_type == :w
                global_dof = LM[1, node_idx]
            elseif dof_type == :θ
                global_dof = LM[2, node_idx]
            end
        else
            # It's the RIGHTMOST node of the very last element
            elem = node_idx - 1
            if dof_type == :w
                global_dof = LM[end-1, elem]
            elseif dof_type == :θ
                global_dof = LM[end, elem]
            end
        end
        
        # C. Apply the dummy equation
        K[global_dof, :] .= 0.0             # 1. Zero out the entire row (note the .= for broadcasting)
        K[global_dof, global_dof] = 1.0     # 2. Set diagonal to 1
        F[global_dof] = val                 # 3. Set Force vector to the prescribed value (usually 0)
    end

    u = K \ F
    
    return (; u, K = K_orig, F = F_orig, mesh, LM)
end

end