function beamParameters()
    beam_length = 8;          # m
    E(x) = 1000^2 * 70;      # GPa -> KN/m^2
    I(x) = 3.5*1e-4;        # m^4
    nElem = -1;
    nnpe = 2;
    return (; beam_length, E, I, nElem, nnpe)
end

function boundaryConditions()
    nBC = 3;
    bcLoc = [0, 0, 4];       # x-location of each support
    bcDOF = [:w, :θ, :w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
    bcVal = [0, 0, 0];       # prescribed value at that DOF

    return (; nBC, bcLoc, bcDOF, bcVal)
end

function forceAndMoments()
    q(x) = if x < 4 0 else 9 end;                # kN/m
    nPointForce = 0;
    nPointMoment = 0;
    pfLoc = [];              # location(s) of point force(s)
    pfVal = [];            # magnitude(s) of point force(s), kN
    pmLoc = [];
    pmVal = [];
    return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
end

function releases()
    nReleases = 0;
    relLoc  = [];       # x-location of release
    relType = [];       # :v (shear release) or :m (moment release)

    return (; nReleases, relLoc, relType)
end
