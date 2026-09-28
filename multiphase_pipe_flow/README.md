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
cell centres, all four field velocities at cell faces. The mesh is a plain
1D non-uniform grid (face positions stored explicitly, no hanging nodes --
see `FlowState.hpp`), starting uniform and optionally locally refined and
coarsened at runtime; see "Adaptive mesh refinement" below. Per time step:

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

### Adaptive mesh refinement

Off by default (`SolverOptions::amr.enabled = false`); when turned on, the
solver periodically (`adaptEveryNSteps`, default every 20 steps) refines
and coarsens the mesh using the **same Kelvin-Helmholtz stability
parameter F already computed for the Andreussi & Persen (1987) interfacial
friction closure** (`Closures.hpp`, `kelvinHelmholtzParameterF()`) as its
indicator: cells where the interface is going unstable (F approaching or
past the closure's own F0 = 0.36 threshold) are split; cells deep in a
stable, quiescent state (F well below threshold) are merged with a like
neighbour. This is the same indicator used by Gourma, Jia & Thompson
(2013), *"Two-Fluid Model for 1D Gas-Liquid Slug Flows: Realizable Mean
Slug Characteristics,"* Multiphase Science and Technology 25(1), 57-79,
for AMR on this class of 1D two-fluid model.

Being a 1D problem, this is a plain non-uniform grid (each cell has
exactly one left and one right neighbour, no hanging nodes or tree
bookkeeping) rather than the block-structured/tree AMR that a multi-D
solver needs -- refining is "insert a face and duplicate the cell's state
into the two halves" and coarsening is "remove a shared face and
width-weight-average the pair"; both are exactly mass-conservative for the
volume fractions. See `FourFieldSolver.cpp`, `adaptMesh()` and
`computeIndicator()`.

By default the mesh is never refined finer than the run's initial spacing
(`minCellWidthFraction = 1.0`) -- AMR's job out of the box is purely to
coarsen quiescent regions relative to a baseline resolution the user
already chose to be adequate (e.g. ~1 diameter per cell, per the paper's
own convergence finding), so it can only reduce cost relative to a fixed
run at that resolution, never increase it by over-refining. Lower
`minCellWidthFraction` explicitly to resolve hot zones finer than the
baseline. On the horizontal slug-formation demo case (`mfs_demo
horizontal_amr`), this gives roughly a **7x speedup** (~300 cells /
~8,200 steps fixed vs. ~65-105 cells / ~3,500 steps adaptive) while
tracking the same interfacial instability (peak F within the run-to-run
noise of the fixed-grid case) -- reproduce with:

```sh
./build/mfs_demo horizontal_amr
```

This h-refinement approach is one of two AMR strategies in the two-phase
pipe-flow literature for this class of model; the other -- a *moving-grid*
scheme, where grid points move with fronts instead of a fixed grid being
locally refined -- is implemented too, see below.

### Moving-mesh (r-adaptive) tracking

Off by default (`SolverOptions::movingMesh.enabled = false`), and
mutually exclusive with `amr` by convention (both change the mesh; combining
them is untested). Where h-refinement keeps a fixed background resolution
and adds/removes cells within a budget, this keeps the cell count **fixed**
and instead continuously relocates every node toward wherever a monitor
function is largest, via the classical equidistribution principle (de Boor,
1974): the mesh is chosen so that the integral of the monitor function is
equal over every cell. A monitor built from the local liquid-holdup
gradient (the classical "arc-length" monitor, Huang & Russell) does the
tracking itself, once a front exists; a smaller term from the same
Kelvin-Helmholtz `F` used by `amr` gives it anticipatory pull toward a
forming front before its holdup gradient is yet sharp. Relocation happens
every step by default (`relocateEveryNSteps`), under-relaxed
(`relaxation`) to avoid the mesh itself oscillating near a still-forming
front, with the new mesh's fields obtained by exact conservative remapping
(cell-centred fields) or linear interpolation (the point-sampled
velocities) from the old one -- see `FourFieldSolver.cpp`,
`computeMonitorFunction()` and `relocateMesh()`.

Because resolution here is limited only by the total point budget N, not
by any minimum-cell-width floor (h-refinement's `minCellWidthFraction`
explicitly forbids refining past the initial spacing by default), this can
concentrate resolution on a sharp front far more tightly than h-refinement
can. On the horizontal slug-formation demo (`mfs_demo
horizontal_movingmesh`), the steepest captured liquid-holdup gradient at
the same simulated time is **~8.7x sharper** than the fixed grid achieves
(down to cells ~0.1 diameter wide right at the front, an order of
magnitude finer than either the fixed grid or AMR's baseline-floor
resolution there). This is **not** a speed optimisation the way AMR is --
relocating and remapping the whole field every step costs more per step
than the fixed grid, and the much finer local cells it creates also force
a smaller Courant-limited time step (more, smaller steps overall); the
same demo case runs roughly 24x *slower* wall-clock for that 8.7x sharper
front. The trade being made is accuracy (less numerically diffused fronts)
for cost, which is exactly the trade the moving-grid/tracking literature
for this model class describes (see Nydal & Banerjee, 1996; De Leebeeck,
2010, "A roll wave and slug tracking scheme for gas-liquid pipe flow").

This implementation is a **periodic-equidistribution r-adaptive** scheme
(relocate the existing point budget every step via de Boor's algorithm),
which is simpler to implement robustly than -- and should not be confused
with -- the object-oriented **Lagrangian slug-tracking** schemes in that
same literature (Nydal & Banerjee 1996; Renault 2007's LASSI), where
individual slugs and bubbles are tracked as discrete objects with their
own integral mass/momentum balances rather than as a locally dense region
of an otherwise-Eulerian field. That alternative is a materially different
model, not just a different mesh strategy, and is not implemented here.

We used this scheme to follow up directly on a finding from
[VALIDATION.md](VALIDATION.md): that the model's interfacial closure
reports several annular-labelled test conditions as unstable (`F` past its
own threshold) while the fixed-grid solution shows no growing wave even at
the paper's own stated converged resolution. Moving-mesh tracking **does**
sharpen a front that is already forming (the 8.7x result above), but on
the specific stuck case from VALIDATION.md it did **not** produce genuine
wave growth from an initially smooth interface -- there was no existing
gradient for the monitor function to concentrate resolution around. That
narrows the likely root cause: it points away from *mesh resolution alone*
and toward the explicit/semi-implicit momentum scheme's treatment of the
disturbance-growth mechanism itself, which is now the more specific target
for further investigation (see VALIDATION.md's "Recommended follow-up").

This is a research/reference-grade implementation: faithful to the
described physics and defensible in every numerical choice it had to make
beyond the paper's text, but it has not been tuned or validated against the
paper's own published figures (that would require the exact closure
coefficients and grid/time-step choices used in the original, proprietary
code). Treat it as a solid, documented starting point for further
calibration, not a certified reproduction of the paper's results.

### Higher-order / flux-limited advection

Off by default (`SolverOptions::advectionLimiter = FluxLimiterType::None`,
i.e. plain first-order upwind, so existing behaviour is unchanged unless
opted into). When set to `Minmod`, `VanLeer`, `Superbee`, or `MC`, two
places switch to a MUSCL/TVD-limited reconstruction instead of plain
upwinding, both built on the same Sweby (1984) "high resolution" form
`recon = donor + 0.5*psi(r)*(acceptor - donor)` (`r` = ratio of the
upwind-side to the local gradient, `psi` the chosen limiter function;
`FluxLimiter.hpp` and the `limitedFaceValue()` helper in
`FourFieldSolver.cpp`):

- The field continuity equations' (`updateContinuity()`) advective face
  fluxes reconstruct `el/ed/eg/eb` at each face directly with this form.
- The layer momentum equations' (`updateLayerMomentum()`) own
  self-advection term `u*du/dz` is *not* already a flux divergence (it is
  non-conservative, "advection form"), so it can't be limited the same
  way directly. Instead it is first recast into the equivalent
  conservative flux form `d(u^2/2)/dz` -- the same identity that relates
  Burgers' equation's conservative and non-conservative forms for smooth
  `u` -- letting `limitedFaceValue()` reconstruct `u` itself at each CELL
  CENTRE (the natural "flux point" for a field stored one index staggered
  at cell faces, exactly like `u1`/`u2` are here), and differencing
  `0.5*uHat^2` between neighbouring cell centres. See
  `buildAdvectiveFlux()` in `FourFieldSolver.cpp` for the full derivation
  in comments. Both are formally second-order accurate in smooth regions
  while remaining Total-Variation-Diminishing (no new overshoot/
  oscillation at a front) -- the limiter falls back to `psi=0` (first
  order) at a local extremum, which is exactly what keeps it bounded.

The pressure-correction system (`solvePressureCorrection()`) deliberately
keeps plain upwind regardless of this setting -- a standard "deferred
correction" split: that step only needs a stable linearization to drive
the pressure iteration, not final transport accuracy, and its coefficients
were tuned against the robustness work above.

First-order upwind is highly numerically diffusive, which is exactly the
mechanism [VALIDATION.md](VALIDATION.md) flags as a candidate explanation
for why some Kelvin-Helmholtz-unstable (F > F0) conditions fail to grow
into the slug/roll-wave regime in this explicit scheme. On the horizontal
slug-formation demo (`mfs_demo horizontal_limiter`), on the *same* fixed
N=300 grid throughout (no mesh change), the steepest captured
liquid-holdup gradient is **~1.8-2.2x sharper** than plain upwind across
the three limiters, at a 1.7-2.4x step-count cost (a steeper resolved
velocity field tightens the Courant-limited time step somewhat):

```sh
./build/mfs_demo horizontal_limiter
```

We also re-ran the specific "stuck" case from VALIDATION.md (D=51mm,
Vsl=0.0025 m/s, Vsg=25 m/s, 1 deg incline) with this enabled -- first with
just the continuity-equation limiting, then with the momentum-equation
extension above added too. Both are informative but neither is a fix. With
either, the back-half liquid-holdup profile develops substantially more
structure than plain upwind ever shows (variation across the back half
goes from ~1-2% to several-hundred percent of the local mean; adding the
momentum-equation term pushes the peak captured gradient roughly another
5x higher again). Tracked over time, though, that structure is the solver
draining down from its initial condition (eL0 = 0.05) toward the
much-thinner film the actual inlet rates imply, resolved more sharply than
upwind's heavy diffusion allowed -- **not** an unboundedly amplifying
wave: the peak gradient rises sharply during that transient, then
saturates (within roughly a factor of 2, in both variants) rather than
continuing to grow once the transient has passed. Both AMR and moving-mesh
(above) already ruled out mesh resolution as the explanation for this
case; this now rules out numerical diffusion in *both* the mass-transport
and the layer-momentum advection schemes too. That leaves the remaining
untried candidates from VALIDATION.md's "Recommended follow-up": the
semi-implicit pressure-velocity coupling's own damping effect on
disturbance growth, and directly measuring the seeded disturbance's actual
growth rate against the inviscid Kelvin-Helmholtz prediction the closure
itself is based on (VALIDATION.md item 1) -- which would also clarify
whether the seeded disturbance's amplitude/frequency is simply a poor
match to this condition's most-unstable wavelength, independent of the
numerical scheme entirely.

### Pressure-velocity coupling investigation

Follow-up on the last remaining numerical-scheme candidate above: could
the segregated, once-per-step pressure correction (`solvePressureCorrection()`,
under-relaxed by `pressureRelaxation`) or the lagged (start-of-step)
interfacial-shear coefficient it and `updateLayerMomentum()` both use
(computed once per step by `computeClosures()`, frozen for that whole
step) be numerically damping the stuck case's disturbance growth? Two
direct experiments, no code changes -- both purely diagnostic, run against
the unmodified solver at its default (first-order) advection scheme:

1. **Time-step refinement.** Any step-size-dependent numerical lag from
   either mechanism should shrink as dt -> 0, recovering more of the true
   growth rate. Instead, the stuck case's peak downstream gradient
   **converges cleanly as dt shrinks**: `courantTarget` 0.5 -> 0.1 -> 0.02
   (a 25x finer step) gives peak gradient 0.01440 -> 0.01427 -> 0.01426
   /m, agreeing to under 0.1%. The time-marching scheme is a *consistent*
   approximation of whatever the discretized system actually does, not a
   source of growth-suppressing lag that fades with resolution.
2. **Direct relaxation sweep.** Varying `pressureRelaxation` at a fixed,
   baseline time step: removing it entirely (0.2 -> 1.0, a full unrelaxed
   correction every step) raises the peak gradient only ~19% (0.0153 ->
   0.0183 /m) -- small, bounded, and nothing like the order-of-magnitude,
   *sustained* amplification that suppressed KH growth would look like.
   `maxPressureChangeFraction` (the per-step correction cap) had no
   measurable effect at any value tested (0.05, 0.2, 1.0 all identical),
   meaning it was never actually binding at this condition.

**This rules out the pressure-velocity coupling as the cause too.** With
mesh resolution, both advection schemes, and now this all ruled out by
direct experiment -- and the solution converging cleanly under both
spatial and temporal refinement -- the likely explanation has shifted away
from "a fixable numerical artifact is suppressing real growth" and toward
this being what these closures' *converged* solution actually does at
this specific condition. The Andreussi & Persen `F` correlation is an
empirical friction-factor enhancement fit to steady wavy-stratified data,
not a first-principles linear-stability growth rate, so `F > F0` flagging
"enhanced interfacial friction" does not necessarily imply that this exact
closure combination's converged solution amplifies a small perturbation
exponentially -- or the seeded disturbance's amplitude/frequency may
simply be a poor match to whatever this system's actual most-unstable mode
is. Directly measuring the disturbance's growth rate against the inviscid
KH prediction (VALIDATION.md item 1) is the most direct way left to tell
those two apart, and is the concrete next step.

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
./build/mfs_demo horizontal             # stratified -> slug transient, cf. Section 4/5.2
./build/mfs_demo vertical               # vertical bubbly riser (extrapolated closures)
./build/mfs_demo terrain                # terrain-following V-section pipeline
./build/mfs_demo horizontal_amr         # fixed-grid vs. h-refinement AMR: speed comparison
./build/mfs_demo horizontal_movingmesh  # fixed-grid vs. moving mesh: front-sharpness comparison
./build/mfs_demo horizontal_limiter     # first-order upwind vs. MUSCL/TVD flux limiters
./build/mfs_demo all                    # runs all six (the mesh-adaptivity ones take a while)
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
