# multiphase_pipe_flow

A C++17 implementation of the **four-field, flow-regime-independent
multi-field model** for transient gas-liquid pipe flow described in:

> M. Bonizzi, P. Andreussi, S. Banerjee, "Flow regime independent, high
> resolution multi-field modelling of near-horizontal gas-liquid flows in
> pipelines", *International Journal of Multiphase Flow* 35 (2009) 34-46.

and **generalised from the paper's near-horizontal scope to arbitrary pipe
inclination**, from horizontal through inclined to fully vertical, including
pipelines whose inclination varies along their length (terrain-following
lines).

## The physical model

Four fields are tracked along the pipe axis: continuous liquid (`el`),
dispersed liquid / droplets (`ed`), continuous gas (`eg`), dispersed gas /
bubbles (`eb`), with `el+ed+eg+eb = 1`. No flow-regime map or flow-regime
dependent closure is ever consulted by the solver; stratified, slug,
annular and bubbly behaviour all emerge from the same equations and the
same closures, exactly as in the source paper. A flow-regime *classifier*
is provided purely as a diagnostic / validation tool (Section 5.1 of the
paper) and has zero influence on the calculation.

Implemented equations (paper equation numbers in parentheses):

* Layer-1 / layer-2 mixture momentum (6)-(7), with the slip-velocity flux
  term dropped, as justified in the paper's Appendix A (shown there to be
  an order of magnitude or more smaller than the terms retained).
* Dispersed-phase (droplet, bubble) momentum (10), reduced to an algebraic
  drag/pressure/gravity balance because the inertial terms were shown to be
  negligible -- exactly the simplification the paper itself makes, and the
  reason the overall system stays hyperbolic at fine grid resolution.
* The four field continuity equations (11)-(14) and the total liquid/gas
  continuity equations (15)-(16).
* The pressure equation (22), obtained by combining the total liquid and
  gas continuity equations and dividing by their respective densities
  (isothermal ideal-gas EOS for the gas, incompressible liquid), solved
  with a standard segregated pressure-correction (SIMPLE-family) scheme,
  per the paper's own citation of Ferziger & Peric (1999).
* All of Tables 1-3's closure relationships: wall and interfacial friction
  factors (Taitel & Dukler 1976 as the default, plus the Table 2
  alternatives -- Spedding & Hand, Moody rough-wall, Andreussi & Persen,
  Andritsos & Hanratty, Cohen & Hanratty, Wallis), bubble
  entrainment/disengagement, droplet entrainment/deposition, and the
  Tomiyama (bubble) / Alipchenkov (droplet) drag coefficients with the
  Andreussi (bubble) / Sarkhi & Hanratty (droplet) diameter closures.
* Eq. (23), the Courant-limited explicit time step.
* Eq. (24)-style mass-conservation bookkeeping (`FourFieldSolver::totalLiquidMass()`,
  `totalGasMass()`, `boundaryMassFluxes()`).

### Generalisation to inclined and vertical pipes

The paper's own momentum equations already carry a general inclination
angle `theta`: a `sin(theta)` gravity-along-the-pipe term and a
`cos(theta)` "hydraulic head" / level-gradient term. This implementation
keeps `theta` as a per-cell field (`FourFieldSolver::setInclinationProfile`),
so the same equations apply unchanged from horizontal (`theta=0`) through
any incline up to vertical (`theta=+-90 deg`), and can vary along the pipe
to represent terrain -- exactly the motivating scenario described in the
paper's own introduction ("flow regimes may change in regions where
pipeline inclination changes due to the terrain").

At `theta = 90 deg` the `cos(theta)` level-gradient term vanishes
identically: there is no longer a preferred low side of the pipe to drive
gravity-segregated stratification, which is the physically correct limit.
**The closures in Tables 1-3, however, were derived and validated by their
original authors for near-horizontal stratified/slug/bubbly/annular flows.**
Applying them unchanged at steep or vertical inclination is an
extrapolation beyond the paper's validated range -- the paper's own
conclusion explicitly flags this extension ("whether a similar approach
might be useful for flows in vertical pipelines") as future work. The
pipe-cross-section geometry (wetted perimeters, interfacial width) is
likewise computed from the classical Taitel & Dukler circular-segment
formulas, which lose their literal physical meaning (there is no
gravity-segregated "bottom" of a vertical pipe) but are retained as an
effective area-based closure, consistent with how most industrial
transient multiphase pipe-flow codes extend a near-horizontal closure set
across the full inclination range. See `include/mfs/PipeGeometry.hpp` for
the detailed discussion. A `vertical` example case is included precisely
to make this extrapolation, and its limits, explicit and runnable.

## Where this implementation had to fill gaps

The source PDF's text extraction lost some minus signs and square-root
radicals (a known artifact of this particular typeset). Two consequences,
both documented in code comments at the point of use:

* The **wall/interfacial friction and gravity source-term signs** in the
  layer momentum equations were reconstructed from first-principles
  physical reasoning (friction opposes motion; interfacial drag
  accelerates the slower phase and decelerates the faster one) rather than
  trusted verbatim from the extracted OCR text, which showed inconsistent
  "+" signs in several places that would be non-physical if taken literally.
* The **droplet entrainment rate** (Pan & Hanratty, 2002) and the
  **interfacial "wave celerity"** feeding the bubble entrainment closure
  could not be recovered as exact closed-form expressions from the OCR
  text. Literature-consistent, clearly-flagged substitutes are used (see
  `Closures.hpp`) and are trivial to replace, since every closure is a
  free function of only physical arguments (the paper's own stated design
  goal: "the software has been formulated in a manner such that closure
  relationships can be readily changed").

## Numerical method

Staggered grid (Harlow & Welch, 1965): volume fractions and pressure at
cell centres, all four field velocities at cell faces. Per time step:

1. Geometry and gas density update from the current pressure/holdup field.
2. Closures (friction factors, shear stresses, entrainment/disengagement
   rates) from the previous step's velocities (lagged coefficients).
3. Layer-1/layer-2 momentum predictor (explicit).
4. Dispersed-phase (droplet, bubble) algebraic momentum solve.
5. Pressure-correction: a tridiagonal (Thomas-algorithm) solve for the
   pressure field that drives the combined liquid+gas mass-flux divergence
   toward zero, with the isothermal ideal-gas compressibility folded in as
   a local accumulation term; velocities are corrected accordingly.
6. Explicit update of the four continuity equations using the final,
   pressure-consistent velocities, with upwind (donor-cell) advection.

Time step is Courant-limited per Eq. (23) (`FourFieldSolver::stableTimeStep()`).

### Numerical robustness safeguards

A few cells getting numerically stiff (not the bulk flow) is where this
kind of solver actually breaks in practice, so several defensive measures
are built in, all in `SolverOptions` and all documented at their point of
use in `FourFieldSolver.cpp`/`.hpp`:

* **Momentum-fraction floor** (`momentumFractionFloor`): the layer-1/2
  momentum equations, and the back-substitution from mixture to
  continuous-phase velocity (Eq. 3-4 inverted), divide by a layer's volume
  fraction. As a layer is pinched toward zero -- e.g. the thin gas film
  just before a liquid slug bridges the pipe, or a vertical riser fed with
  almost no continuous gas -- that division can amplify an otherwise
  ordinary shear or entrainment term far past what's physical. A dedicated,
  coarser floor is applied only to these denominators (not to the tracked
  volume fraction itself).
* **Velocity limiter** (`maxVelocity`): a generous but finite bound on all
  four field velocities, well above realistic pipe-flow speeds, absorbing
  transient overshoots from the explicit/lagged treatment of locally stiff
  drag and entrainment terms before they can compound.
* **Pressure-correction relaxation and clamp** (`pressureRelaxation`,
  `maxPressureChangeFraction`, `minPressure`/`maxPressureFactor`): the
  segregated pressure solve is a single linear pass per time step, not an
  iterated SIMPLE loop, so its correction is under-relaxed and clamped
  (both a per-step fractional cap and an absolute floor/ceiling) rather
  than applied in full, which is what a cell with a very weak local
  compressibility coefficient (near-zero gas fraction) needs to stay
  bounded.
* The dispersed-phase (droplet/bubble) algebraic momentum balance (Eq. 10)
  is deliberately re-derived in `updateDispersedMomentum()` (see the
  comment there) so that its denominator is the **entrainment** rate, not
  the disengagement/deposition rate -- the latter cancels out of the
  velocity equation by the quotient rule, and using it instead produces a
  spurious 0/0 whenever a dispersed fraction is locally zero but actively
  being fed by entrainment.

None of these change the equations being solved in the bulk of the flow;
they only bound the response in the small number of cells where the
discretisation is locally stiff, and are conservative enough that the
horizontal demo case's slug-formation transient is essentially unaffected
by them (verified by running with tighter/looser settings and comparing).

This is a research/reference-grade implementation: faithful to the
described physics and defensible in every numerical choice it had to make
beyond the paper's text, but it has not been tuned or validated against the
paper's own published figures (that would require the exact closure
coefficients and grid/time-step choices used in the original, proprietary
code). Treat it as a solid, documented starting point for further
calibration, not a certified reproduction of the paper's results.

## Validation against experimental data

**[VALIDATION.md](VALIDATION.md)** compares the solver's predicted flow
regime against Shoham's (1982) classic experimental dataset spanning the
full -90 deg to +90 deg inclination range (182 cases, zero numerical
failures). Headline result: **91.7% agreement on stratified-labeled
conditions across every inclination tested (100% horizontal, 85%
inclined)**, but poor agreement (0%) on annular/slug conditions, traced to
a specific, reproducible finding -- the interfacial friction closure
itself reports these conditions as unstable, but the explicit time-
marching solution isn't amplifying the seeded disturbance into a growing
wave. Read VALIDATION.md for the full methodology, breakdown, and
recommended follow-up. The validation driver (`validate_shoham`, built by
default) and its input/output CSVs are in `validation/`.

## Building

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Requires a C++17 compiler and CMake >= 3.15. No external dependencies.

## Running

```sh
./build/mfs_demo horizontal   # stratified -> slug transient, cf. Section 4/5.2
./build/mfs_demo vertical     # vertical bubbly riser (extrapolated closures)
./build/mfs_demo terrain      # terrain-following V-section pipeline
./build/mfs_demo all          # runs all three
```

Each case writes a CSV time series (`<case>_output.csv`) of cell-centred
fields (`eL, eG, el, ed, eg, eb, P, u1, u2, ul, ug, ud, ub, theta`) and
prints a flow-regime-cell-count summary and a mass-balance diagnostic to
stdout at each snapshot.

## Library usage

```cpp
#include "mfs/FourFieldSolver.hpp"

mfs::FluidProperties fluid;               // defaults: air/water at ~20 C
mfs::FourFieldSolver solver(/*D=*/0.08, /*L=*/30.0, /*N=*/300, fluid);
solver.setInclinationConstant(0.0);       // horizontal; or any angle in radians

mfs::BoundaryConditions bc;
bc.inletSuperficialLiquid = 0.2;          // [m/s]
bc.inletSuperficialGas = 3.0;             // [m/s]
bc.inletLiquidHoldup = 0.15;
bc.outletPressure = 1.0e5;                // [Pa]
solver.setBoundaryConditions(bc);
solver.initializeStratified(0.15);

while (solver.time() < 10.0) {
    solver.step(solver.stableTimeStep());
}
```

## Directory layout

```
include/mfs/    Public headers (one class/module per header)
src/            Implementation + CLI driver (main.cpp)
```

| File | Contents |
|---|---|
| `PipeGeometry.hpp/.cpp` | Circular-pipe stratified-layer geometry |
| `FluidProperties.hpp` | Liquid (incompressible) / gas (ideal-gas) properties |
| `Closures.hpp/.cpp` | Tables 1-3: friction, entrainment, drag, particle size |
| `FlowState.hpp` | Staggered-grid field storage |
| `TridiagonalSolver.hpp` | Thomas algorithm |
| `FourFieldSolver.hpp/.cpp` | The solver itself |
| `FlowRegimeClassifier.hpp/.cpp` | Diagnostic-only Section 5.1 regime criteria |
| `main.cpp` | CLI driver / example cases |
