module preprocessor

export generateMesh,generateLM;

function generateMesh(beam,bc, fm, rel;)
    specialLocs = vcat(bc.bcLoc, fm.pfLoc, fm.pmLoc, rel.relLoc)
    specialLocs = sort(unique(specialLocs))

    minElements = length(specialLocs) - 1
    nElem = max(beam.nElem, minElements)
    
    if beam.nElem < minElements
        @warn "Requested $(beam.nElem) elements; bumping to $minElements to guarantee a node at every support/force/moment/release."
    end

    extra = nElem - minElements                     
    base, remainder = divrem(extra, minElements)
    subdivPerSeg = fill(base + 1, minElements)       
    subdivPerSeg[1:remainder] .+= 1                  

    nodeLocs = Float64[specialLocs[1]]
    for s in 1:minElements
        pts = range(specialLocs[s], specialLocs[s+1], length = subdivPerSeg[s] + 1)
        append!(nodeLocs, pts[2:end])               
    end

    releaseElemIdx  = [findfirst(==(loc), nodeLocs) for loc in rel.relLoc]
    releaseNodeType = rel.relType

    return (; nodeLocs, nElem, releaseElemIdx, releaseNodeType)
end

function generateLM(mesh)
    nElem = mesh.nElem
    releaseElemIdx = mesh.releaseElemIdx
    releaseNodeType = mesh.releaseNodeType

    LM = zeros(Int, 4, nElem)
    current_dof = 1

    LM[1, 1] = current_dof; current_dof += 1
    LM[2, 1] = current_dof; current_dof += 1
    LM[3, 1] = current_dof; current_dof += 1
    LM[4, 1] = current_dof; current_dof += 1

    for i in 2:nElem
        # Check if the current element's left node has a release
        rel_idx = findfirst(==(i), releaseElemIdx)

        if isnothing(rel_idx)
            # NO RELEASE: Stitch both left DOFs to the previous element's right DOFs
            LM[1, i] = LM[3, i-1]
            LM[2, i] = LM[4, i-1]
        else
            # RELEASE FOUND: Split one DOF, stitch the other
            rel_type = releaseNodeType[rel_idx]
            
            if rel_type == :m
                # Moment release: share deflection (w), split rotation (theta)
                LM[1, i] = LM[3, i-1]                      # Shared shear
                LM[2, i] = current_dof; current_dof += 1   # Independent moment DOF
            elseif rel_type == :v
                # Shear release: split deflection (w), share rotation (theta)
                LM[1, i] = current_dof; current_dof += 1   # Independent shear DOF
                LM[2, i] = LM[4, i-1]                      # Shared moment
            end
        end

        # The right node of the current element always gets new DOFs
        LM[3, i] = current_dof; current_dof += 1
        LM[4, i] = current_dof; current_dof += 1
    end

    return LM
end


end