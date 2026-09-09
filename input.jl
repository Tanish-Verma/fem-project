# our Group
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

# Group 1
# function beamParameters()
#     beam_length = 4;          # in m
#     E(x) = 20* 10^3;      # in GPa -> KN/m^2
#     I(x) = 1;        # in m^4
#     nElem = 4; # number of elements
#     nnpe = 2; # number of nodes per element
#     return (; beam_length, E, I, nElem, nnpe)
# end

# function boundaryConditions()
#     nBC = 5; # number of boundary conditions
#     bcLoc = [0 ,0, 2, 4,4];       # x-location of each support
#     bcDOF = [:w, :θ, :w,:θ ,:w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
#     bcVal = [0, 0, 0,0,0];       # prescribed value at that DOF

#     return (; nBC, bcLoc, bcDOF, bcVal)
# end

# function forceAndMoments()
#     q(x) = if x >= 3 24 else 0 end;                # in kN/m
#     nPointForce = 0;          # number of point forces
#     nPointMoment = 0;         # number of point moments
#     pfLoc = [];              # location(s) of point force(s)
#     pfVal = [];            # magnitude(s) of point force(s), in kN
#     pmLoc = [];              # location(s) of point moment(s)
#     pmVal = [];              # magnitude(s) of point moment(s), in kN m
#     return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
# end

# function releases()
#     nReleases = 1;         # number of releases
#     relLoc  = [3];          # x-location of release
#     relType = [:m];          # :v (shear release) or :m (moment release)

#     return (; nReleases, relLoc, relType)
# end

# Group 9
# function beamParameters()
#     beam_length = 20;          # in m
#     E(x) = 210* 10^6;      # in GPa -> KN/m^2
#     I(x) = (if x < 10 54 else 18 end)* 1e-5;        # in m^4
#     nElem = 4; # number of elements
#     nnpe = 2; # number of nodes per element
#     return (; beam_length, E, I, nElem, nnpe)
# end

# function boundaryConditions()
#     nBC = 3; # number of boundary conditions
#     bcLoc = [0,10,20];       # x-location of each support
#     bcDOF = [:w,:w,:w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
#     bcVal = [0, 0, 0];       # prescribed value at that DOF

#     return (; nBC, bcLoc, bcDOF, bcVal)
# end

# function forceAndMoments()
#     q(x) = if x >= 10 36 else 0 end;                # in kN/m
#     nPointForce = 1;          # number of point forces
#     nPointMoment = 0;         # number of point moments
#     pfLoc = [5];              # location(s) of point force(s)
#     pfVal = [360];            # magnitude(s) of point force(s), in kN
#     pmLoc = [];              # location(s) of point moment(s)
#     pmVal = [];              # magnitude(s) of point moment(s), in kN m
#     return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
# end

# function releases()
#     nReleases = 0;         # number of releases
#     relLoc  = [];          # x-location of release
#     relType = [];          # :v (shear release) or :m (moment release)

#     return (; nReleases, relLoc, relType)
# end

#random grp
# function beamParameters()
#     beam_length = 15;          # in m
#     E(x) = 210* 10^6;      # in GPa -> KN/m^2
#     I(x) = 600* 1e-6;        # in m^4
#     nElem = 4; # number of elements
#     nnpe = 2; # number of nodes per element
#     return (; beam_length, E, I, nElem, nnpe)
# end

# function boundaryConditions()
#     nBC = 3; # number of boundary conditions
#     bcLoc = [0,0,10];       # x-location of each support
#     bcDOF = [:w,:θ,:w];    # which DOF is restrained: :w (deflection) or :θ (rotation)
#     bcVal = [0, 0, 0];       # prescribed value at that DOF

#     return (; nBC, bcLoc, bcDOF, bcVal)
# end

# function forceAndMoments()
#     q(x) = if x >= 10 12 else 0 end;                # in kN/m
#     nPointForce = 0;          # number of point forces
#     nPointMoment = 0;         # number of point moments
#     pfLoc = [];              # location(s) of point force(s)
#     pfVal = [];            # magnitude(s) of point force(s), in kN
#     pmLoc = [];              # location(s) of point moment(s)
#     pmVal = [];              # magnitude(s) of point moment(s), in kN m
#     return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
# end

# function releases()
#     nReleases = 0;         # number of releases
#     relLoc  = [];          # x-location of release
#     relType = [];          # :v (shear release) or :m (moment release)

#     return (; nReleases, relLoc, relType)
# end