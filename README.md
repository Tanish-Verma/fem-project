# FEM Beam Solver

This repository contains a Julia implementation of the finite element method
(FEM) for Euler-Bernoulli beams. It supports:

- spatially varying Young's modulus `E(x)` and second moment of area `I(x)`;
- distributed loads, point forces, and point moments;
- prescribed displacement (`:w`) and rotation (`:θ`) boundary conditions;
- moment releases (`:m`) and shear releases (`:v`);
- two-node and higher-order beam elements; and
- displacement, rotation, shear-force, and bending-moment post-processing.

The default example in `input.jl` uses metres, kN, kN/m, kN m, and kN/m².
Keep all input quantities in one consistent unit system.

## Requirements

- Julia 1.12.x. The included manifest was generated with Julia 1.12.7.
- A terminal or the Julia REPL.
- A graphical environment if you want to display plots interactively. Plot
  files are also saved to disk, so interactive display is not required for
  convergence studies.

The repository is a Julia project. `Project.toml` describes the direct
dependencies and `Manifest.toml` records the resolved dependency versions.

## Download and install

Clone the repository and enter its directory:

```bash
git clone https://github.com/Tanish-Verma/fem-project.git
cd fem-project
```

Instantiate the project environment. This downloads all packages recorded by
the project and manifest:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

You can verify that Julia sees the project with:

```bash
julia --project=. -e 'using Pkg; Pkg.status()'
```

On a fresh machine, the first run can take longer because Julia compiles
packages and the plotting backend.

## Run the beam example

From the repository root, run:

```bash
julia --project=. Beam.jl
```

`Beam.jl` includes `input.jl` and the `FEMBeamSolver` module, then performs
these operations:

1. Builds the beam, boundary-condition, load, and release data from
   `input.jl`.
2. Calls `Solve`, which generates the mesh and location matrix, assembles the
   global stiffness matrix and load vector, applies boundary conditions, solves
   for the displacement vector, and computes reactions.
3. Prints the displacement vector `result.u` and reaction vector
   `result.reactions`.
4. Calls `postprocess` and `plot_results`.
5. Saves the combined plot as `combined_results.pdf` and displays it. The
   script waits for Enter before closing the plot window.

The displacement vector interleaves deflection and rotation degrees of freedom
according to the generated location matrix. The exact length depends on the
number of elements, polynomial order, and releases.

## How to change a problem?

Edit the active definitions in `input.jl`. Each function returns a named tuple
that is consumed by the solver.

### Beam parameters

```julia
function beamParameters()
  beam_length = 8
  E(x) = 210 * 1e6
  I(x) = 100 * 1e-6
  nElem = 4
  nnpe = 2
  return (; beam_length, E, I, nElem, nnpe)
end
```

- `beam_length`: total beam length.
- `E(x)`: Young's modulus as a function of position.
- `I(x)`: second moment of area as a function of position.
- `nElem`: requested number of elements.
- `nnpe`: nodes per element. `2` is the standard two-node cubic Hermite
  beam element; larger values enable higher-order interpolation.

`E` and `I` must be functions, even when the value is constant. For example,
use `E(x) = 210e6`, not just `E = 210e6`.

### Boundary conditions

```julia
function boundaryConditions()
  nBC = 3
  bcLoc = [0, 4, 8]
  bcDOF = [:w, :w, :w]
  bcVal = [0, 0, 0]
  return (; nBC, bcLoc, bcDOF, bcVal)
end
```

Each entry describes one prescribed degree of freedom:

- `nBC` is the number of Boundary Conditions.
- `bcLoc[i]` is its x-coordinate and must be a beam node location.
- `bcDOF[i]` is `:w` for deflection or `:θ` for rotation.
- `bcVal[i]` is the prescribed value, normally zero.

The mesh always creates nodes at supports, point-load locations, point-moment
locations, and release locations. If `nElem` is too small to do this, the
preprocessor increases it automatically and emits a warning.

### Distributed and concentrated loads

```julia
function forceAndMoments()
  q(x) = 90
  nPointForce = 1
  nPointMoment = 0
  pfLoc = [6]
  pfVal = [120]
  pmLoc = []
  pmVal = []
  return (; q, nPointForce, nPointMoment, pfLoc, pfVal, pmLoc, pmVal)
end
```

- `q(x)` is the distributed load function.
- `nPointForce`and `nPointMoment`are the number of point forces and moments applied on the Beam.
- `pfLoc` and `pfVal` contain matching point-force positions and magnitudes.
- `pmLoc` and `pmVal` contain matching point-moment positions and magnitudes.
- Keep the location/value arrays the same length. The solver uses the paired
  entries with `zip`.

### Releases

```julia
function releases()
  nReleases = 0
  relLoc = []
  relType = []
  return (; nReleases, relLoc, relType)
end
```

For each release, put the location in `relLoc` and use the corresponding entry
in `relType`:

- `:m`: moment release; rotation is split while deflection remains shared.
- `:v`: shear release; deflection is split while rotation remains shared.

The active input file contains several commented example groups. To use one,
replace or uncomment the active four function definitions, rather than leaving
multiple definitions of the same function active at once.

## Convergence study

Run:

```bash
julia --project=. ConvergenceStudy.jl
```

The script uses the active problem in `input.jl` and performs:

- an h-refinement study by increasing the element count while keeping
  polynomial degree fixed; and
- a p-refinement study by increasing `nnpe` while keeping the base mesh.

It saves `question6_convergence.png` and `question7_convergence.png`. In an
interactive terminal it displays each plot and waits for Enter between plots.

## Tests

Run the complete test suite from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover element stiffness and load vectors, higher-order elements,
mesh and location-matrix generation, releases and concentrated loads, the
solver, post-processing, and h/p convergence metrics.

## Project layout

### Top-level files

- `Project.toml`: Julia project metadata, direct dependencies, compatibility
  bounds, and the test target.
- `Manifest.toml`: machine-generated, exact dependency resolution for the
  project. Do not edit it by hand; regenerate it with Julia's package manager
  when dependencies change.
- `input.jl`: the user-editable beam model and loading definition. It defines
  `beamParameters`, `boundaryConditions`, `forceAndMoments`, and `releases`.
- `Beam.jl`: the main runnable example and normal entry point.
- `ConvergenceStudy.jl`: the h- and p-refinement experiment and plot writer.
- `README.md`: this usage and architecture guide.
- `combined_results.pdf`: a saved four-panel result plot for a previously run
  example.
- `question6_convergence.png`: a previously saved h-refinement plot.
- `question7_convergence.png`: a previously saved p-refinement plot.
- `project_problems_2026 page 7.pdf`: an external problem/reference document
  stored with the project; it is not loaded by the Julia program.

### `src/`

The source directory contains the reusable solver implementation. The files
are included in dependency order by `src/FEMBeamSolver.jl`.

- `FEMBeamSolver.jl`: top-level module that assembles and exports the public
  solver API: `Solve`, `postprocess`, `plot_results`, mesh utilities, element
  routines, quadrature, and shape functions.
- `Polynomial.jl`: small polynomial type with arithmetic, evaluation, and
  differentiation support.
- `GaussQuadrature.jl`: Legendre polynomial construction and Gauss-Legendre
  points and weights.
- `ShapeFunct.jl`: constructs beam shape functions and their derivatives on the
  reference element.
- `Kelem.jl`: computes an element stiffness matrix from `E(x)`, `I(x)`, the
  element length, and interpolation order.
- `Felem.jl`: computes an element load vector from the distributed load `q(x)`.
- `Global_KF.jl`: assembles element matrices and vectors into global `K` and
  `F`, including point forces and point moments.
- `preprocessor.jl`: creates a mesh containing all special locations and builds
  the element-to-global degree-of-freedom location matrix `LM`, including
  release connectivity.
- `Solver.jl`: performs the five-step solve: mesh/LM generation, assembly,
  boundary-condition application, linear solve, and reaction calculation.
- `PostProcessing.jl`: samples deflection and rotation, reconstructs shear and
  moment diagrams, and creates/saves the four-panel plot.
- `ConvergenceMetrics.jl`: computes a sampled maximum absolute deflection for
  convergence studies.

### `test/`

- `runtests.jl`: test entry point loaded by `Pkg.test()`.
- `testpreprocessor.jl`: mesh and `LM` tests, including releases and higher
  order elements.
- `testglobalkf.jl`: global assembly tests for point loads, point moments, and
  releases.
- `testSolver.jl`: end-to-end one-element solver test.
- `testpostprocessing.jl`: shear, moment, reaction, and sampled output tests.
- `testconvergence.jl`: convergence metric and h/p-refinement tests.

## Using the solver from another Julia script

From the repository root, a separate script can reuse the public module:

```julia
include("input.jl")
include("src/FEMBeamSolver.jl")
using .FEMBeamSolver

result = Solve(beamParameters(), boundaryConditions(), forceAndMoments(), releases())
processed = postprocess(result, beamParameters(), boundaryConditions(), forceAndMoments())
plots = plot_results(processed)
```

Run that script with `julia --project=. path/to/script.jl`. The top-level module
is included locally, so this repository does not need to be installed as a
registered package before use.
