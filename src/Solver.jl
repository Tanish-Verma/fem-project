module Solver

include("Global_KF.jl")
include("preprocessor.jl")

using .preprocessor
using .Global_KF
using LinearAlgebra

export Solve

function Solve(beam, bc, fm, rel)
    total_start = time()

    println("Step 1: Generating mesh and location matrix")
    step_start = time()
    mesh = generateMesh(beam, bc, fm, rel)
    LM = generateLM(mesh,beam.nnpe)
    println("Step 1 completed in $(time() - step_start) seconds")
    
    println("Step 2: Assembling global K and F")
    step_start = time()
    K, F = global_kf(mesh, LM, beam, fm, nnpe=beam.nnpe)
    K_orig = copy(K);
    F_orig = copy(F);
    println("Step 2 completed in $(time() - step_start) seconds")

    println("Step 3: Applying boundary conditions")
    step_start = time()
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
    println("Step 3 completed in $(time() - step_start) seconds")

    println("Step 4: Solving the system")
    step_start = time()
    u = K \ F
    println("Step 4 completed in $(time() - step_start) seconds")
    
    println("Step 5: Computing reactions")
    step_start = time()
    reactions = K_orig * u - F_orig
	println("Step 5 completed in $(time() - step_start) seconds")
	println("Total solver time: $(time() - total_start) seconds")
	    
    return (; u, K = K_orig, F = F_orig, reactions, mesh, LM)
end

end
