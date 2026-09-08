function beamParameters()
    beam_length = 4;          # m
    E(x) = 20*10^3;      # GPa -> KN/m^2
    I(x) = 1;        # m^4
    nElem = -1;
    nnpe = 2;
    return (; beam_length, E, I, nElem, nnpe)
end

function boundaryConditions()
    nBC = 5;
    bcLoc = [0, 0, 2, 4, 4];       # x-location of each support
    bcDOF = [:w, :θ, :w, :w, :θ];    # which DOF is restrained: :w (deflection) or :θ (rotation)
    bcVal = [0, 0, 0, 0, 0];       # prescribed value at that DOF

    return (; nBC, bcLoc, bcDOF, bcVal)
end

function forceAndMoments()
    q(x) = if x < 3 0 else 24 end;                # kN/m
    nPointForce = 0;
    nPointMoment = 0;
    pfLoc = [];              # location(s) of point force(s)
    pfVal = [];            # magnitude(s) of point force(s), kN
    pmLoc = [];
    pmVal = [];
    return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
end

function releases()
    nReleases = 1;
    relLoc  = [3];       # x-location of release
    relType = [:m];       # :v (shear release) or :m (moment release)

    return (; nReleases, relLoc, relType)
end
