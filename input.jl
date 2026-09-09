function beamParameters()
    beam_length = 8;          # in m
    E(x) = 210 * 1e6;      # in GPa -> KN/m^2
    I(x) = 100 * 1e-6;        # in m^4
    nElem = 4; # number of elements
    nnpe = 2; # number of nodes per element
    return (; beam_length, E, I, nElem, nnpe)
end

function boundaryConditions()
    nBC = 3; # number of boundary conditions
    bcLoc = [0, 4, 8];       # x-location of each support
    bcDOF = [:w, :w, :w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
    bcVal = [0, 0, 0];       # prescribed value at that DOF

    return (; nBC, bcLoc, bcDOF, bcVal)
end

function forceAndMoments()
    q(x) = 90;                # in kN/m
    nPointForce = 1;          # number of point forces
    nPointMoment = 0;         # number of point moments
    pfLoc = [6];              # location(s) of point force(s)
    pfVal = [120];            # magnitude(s) of point force(s), in kN
    pmLoc = [];              # location(s) of point moment(s)
    pmVal = [];              # magnitude(s) of point moment(s), in kN m
    return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
end

function releases()
    nReleases = 0;         # number of releases
    relLoc  = [];          # x-location of release
    relType = [];          # :v (shear release) or :m (moment release)

    return (; nReleases, relLoc, relType)
end
