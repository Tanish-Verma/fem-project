function beamParameters()
    beam_length = 8;          # m
    E(x) = 1000^2 * 210;      # GPa -> KN/m^2
    I(x) = 100 * 1e-6;        # m^4
    nElem = 3;
    return (; beam_length, E, I)
end

function boundaryConditions()
    nBC = 3;
    bcLoc = [0, 4, 8];       # x-location of each support
    bcDOF = [:w, :w, :w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
    bcVal = [0, 0, 0];       # prescribed value at that DOF

    return (; nBC, bcLoc, bcDOF, bcVal)
end

function forceAndMoments()
    q(x) = 90;                # kN/m
    nPointForce = 1;
    nPointMoment = 0;
    pfLoc = [6];              # location(s) of point force(s)
    pfVal = [120];            # magnitude(s) of point force(s), kN
    pmLoc = [];
    pmVal = [];
    return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
end

function releases()
    nReleases = 1;
    relLoc  = [2];       # x-location of release
    relType = [:m];       # :v (shear release) or :m (moment release)

    return (; nReleases, relLoc, relType)
end
