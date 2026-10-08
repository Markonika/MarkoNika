# Validation against experimental data

This document records an independent validation pass of the solver in
this repository against published experimental two-phase pipe flow data,
covering both the paper's own near-horizontal scope and the inclined /
vertical extension. It is in addition to, not a replacement for, the
comparisons already described in the source paper itself (Section 5),
which this codebase does not attempt to reproduce bit-for-bit (see
README.md, "Numerical method").

## Data source

**Shoham (1982)**, PhD thesis, "Flow pattern transition and
characterization in gas-liquid two-phase flow in inclined pipes", Tel
Aviv University. Air-water flow, pipe diameters 25.4 mm and 51 mm,
inclination swept from **-90 deg (vertical downward) through horizontal to
+90 deg (vertical upward)**, with six flow patterns reported: stratified
smooth (SS), stratified wavy (SW), intermittent/slug (I), annular (A),
dispersed bubble (DB), and bubble (B).

The raw data (5,675 points: `Vsl, Vsg, VisL, VisG, DenL, DenG, ST, Ang,
ID, FlowPattern`) was obtained from a public compilation:

> BioAITeam, "Machine learning applications to predict two-phase flow
> patterns", https://github.com/BioAITeam/Machine-learning-applications-to-predict-two-phase-flow-patterns
> (`Databases/ShohamDB.csv`), accompanying Arteaga-Arteaga et al. (2021),
> "Machine learning applications to predict two-phase flow patterns",
> *PeerJ Computer Science* 7:e798, https://doi.org/10.7717/peerj-cs.798

This dataset was chosen specifically because it is (a) the classic,
widely-cited experimental campaign for flow-pattern transitions across the
*entire* inclination range, which is exactly the horizontal-to-vertical
extension this codebase implements, and (b) freely available as a clean
CSV, letting the comparison be automated rather than hand-transcribed
from a paper's figures.

Other candidate sources found but not used here (either not open-data, or
duplicative of what Shoham already covers): Nydal, Pintus & Andreussi
(1992) and Andritsos, Williams & Hanratty (1989) -- both already used by
the *source paper itself* for slug length/frequency and the
stratified-slug transition boundary respectively, so they validate the
same horizontal regime the paper's own Section 5 already covers, published
only as in-paper tables/figures, not an open dataset; Al-Kayiem et al.,
"Experimental data for the slug two-phase flow characteristics in
horizontal pipeline", *Data in Brief* 16 (2018) 121-125,
https://doi.org/10.1016/j.dib.2017.11.026 (horizontal only, slug
length/velocity, would be a good target for a future quantitative
follow-up); and a 2025 synchronized-sensor horizontal-pipe dataset
(stratified/slug/dispersed bubble, UNICAMP/LabPetro), also horizontal
only. If a genuinely inclined *quantitative* (not just regime-label)
dataset is needed later, Barnea, Shoham, Taitel & Dukler, "Gas-liquid flow
in inclined tubes" (several 1980-85 papers) is the next place to look, but
its data exists only as published tables.

## Method

`validation/validate_shoham.cpp` (built alongside the main library) reads
a CSV of `(Vsl, Vsg, fluid properties, Ang, ID, FlowPattern)` rows, and
for each one:

1. Builds a `FourFieldSolver` with the given diameter, fluid properties
   (gas density is reproduced at atmospheric outlet pressure via the
   ideal-gas constant, since the solver tracks pressure/EOS rather than
   taking a fixed gas density directly), and inclination.
2. Runs from a stratified initial condition with a small seeded
   disturbance for several pipe residence times (capped for tractability
   -- see "Limitations" below).
3. Classifies the resulting flow regime with the library's own
   diagnostic-only classifier (`mfs::classifyFromHistory`, Section 5.1
   criteria), using the last third of the pipe (nearest the outlet, most
   developed) and the back half of the time history (letting the initial
   transient settle).
4. Maps both the experimental label and the model's prediction onto a
   common 4-category scheme (`stratified`, `annular`, `slug`, `bubbly`:
   SS/SW -> stratified, I -> slug, A -> annular, DB/B -> bubbly) and
   records agreement.

182 points were sampled (up to 2 per distinct angle/flow-pattern
combination in the dataset) to keep the sweep's runtime reasonable while
still covering all 23 distinct angles present, from -90 to +90 deg
(`validation/shoham_sample.csv`, `validation/shoham_results_fast_sweep.csv`
for the full per-case output). A 24-case subset of the failures was then
re-run with a much longer pipe (250 diameters vs. 60), finer mesh (100
cells vs. 40, i.e. 1 diameter per cell -- the resolution the source paper
itself identifies as needed for grid-independence) and longer run time, to
separate "genuinely wrong" from "just needed more development length /
resolution" (`validation/shoham_results_long_pipe_subset.csv`).

To reproduce:

```sh
cmake --build build --target validate_shoham   # or compile validation/validate_shoham.cpp
                                                # directly against the mfs sources, see below
./build/validate_shoham validation/shoham_sample.csv /tmp/out.csv
# optional long-pipe args: lengthInDiameters cellsPerRun maxResidenceTimes maxSimTime
./build/validate_shoham /tmp/hard_cases.csv /tmp/out2.csv 250 100 20 20
```

## Results

**Zero numerical failures.** All 182 cases (plus the 24-case long-pipe
subset, 206 runs total) completed without producing NaN/Inf anywhere in
the domain, across the full -90 to +90 degree sweep and a wide range of
superficial velocities (0.001-26 m/s). This is itself a meaningful result:
it's a direct stress test of the inclination generalisation and the
numerical robustness safeguards added during development (see README.md,
"Numerical robustness safeguards") on 206 independent, experimentally
grounded operating points, not just the three curated demo cases.

**Overall flow-regime classification accuracy: 50/182 = 27.5%.** Broken
down, this splits sharply by regime:

| Labeled regime | Accuracy |
|---|---|
| Stratified (SS/SW) | **33/36 = 91.7%** |
| Bubbly (DB/B) | 17/54 = 31.5% |
| Annular (A) | 0/46 = 0% |
| Slug (I) | 0/46 = 0% |

and by inclination, restricted to the stratified-labeled cases where the
model performs well:

| Inclination | Stratified accuracy |
|---|---|
| Horizontal (\|angle\| <= 2 deg) | **16/16 = 100%** |
| Inclined (2 < \|angle\| < 85 deg) | 17/20 = 85% |

Overall accuracy by inclination bucket (all regimes) was 42.3% horizontal,
22.2% inclined, 18.2% near-vertical -- the drop is driven entirely by the
annular/slug/bubbly misses above, which occur at every angle, not
specifically at steep incline; see below.

**The long-pipe/fine-mesh re-test did not change the picture**: of 24
annular/slug failures re-run at 250 diameters and 1-diameter mesh spacing
(vs. the fast sweep's 60 diameters / 2.5-diameter spacing), 0/24 flipped
to correct, still overwhelmingly predicted "stratified". This rules out
*mesh under-resolution alone* as the explanation.

**Root-cause finding (original, 2024).** A representative failing case
(D = 51 mm, Vsl = 0.0025 m/s, Vsg = 25 m/s, 1 deg incline, labeled annular)
was inspected directly: even at 1-diameter mesh spacing and 10+ pipe
residence times, the liquid film height along the back half of the pipe is
essentially flat (1-2% variation -- numerical noise, not a growing wave).
Yet querying the solver's own interfacial-friction closure at the same
local conditions gives F = 0.66 against the Andreussi & Persen (1987)
threshold F0 = 0.36 (see `Closures.hpp`), i.e. **the closure itself
reports this point as being in the enhanced-shear, unstable regime** (an
11.9x enhancement of the interfacial friction factor over its smooth-wall
value) -- but the explicit time-marching solution is not amplifying the
seeded disturbance into a growing wave. That combination (unstable-per-
closure, stable-in-practice) points to the disturbance not growing fast
enough under this codebase's explicit/semi-implicit momentum treatment,
rather than to a wrong closure choice. This is a genuine, reproducible
limitation of the current implementation, not a data or methodology
artifact, and is the most concrete lead for follow-up work.

**Resolution (see "Recommended follow-up" item 1, Update 5, below): this
was a mis-framed comparison, not a solver deficiency.** F0 = 0.36 is
Andreussi & Persen's *empirical* threshold for onset of enhanced
interfacial friction (a wavy, rough-but-not-necessarily-exponentially-
growing interface) -- it is not the threshold for actual inviscid
Kelvin-Helmholtz exponential wave growth, which is F = 1 (recoverable from
this model's own governing equations; see below). At the case's inlet
condition (eL = 0.05) F = 1.14, genuinely above the F = 1 growth
threshold, and the solver does grow the disturbance there, fast, exactly
as it should. But Vsl = 0.0025 m/s is small enough that the film's true
equilibrium thickness is much less than the arbitrary eL = 0.05 starting
point, and as it drains toward that equilibrium (the same "draining
transient" seen throughout this file's mesh/advection/pressure-coupling
follow-up work), F **drops below 1 well before reaching the F = 0.66 point
quoted above** (the crossing is at eL ~ 0.027) -- so by the time the film
reaches the thickness this finding actually queried, growth has correctly
stopped, because the flow is no longer inviscid-unstable there, even
though it is still past the *different* F0 threshold. The flat back-half
profile was the model getting this right, not failing to grow a wave it
should have grown.

## Interpretation

- The **flow-regime-independent stratified prediction is strongly
  validated across the entire inclination range this codebase adds**
  (100% on horizontal, 85% on inclined, using the *same* closures at every
  angle, no regime-specific tuning) -- directly supporting the paper's
  central claim and the specific horizontal/inclined comparison asked for.
- **Wave growth into annular/slug/bubbly regimes is under-*predicted by
  the regime classifier's comparison***, at least for the representative
  case this file traced to ground -- but see the Update below: this
  turned out to be a mis-framed comparison rather than a solver defect,
  for that case. The classifier miss (labeling it "stratified" when
  Shoham's data says "annular") is still real and still stands.
- The original picture here was: *the model correctly stays stratified
  when the flow should be stratified, and correctly resolves stratified
  behaviour identically well from horizontal through steep incline, but
  currently under-predicts the growth of disturbances into the other
  regimes.* Follow-up work (below) revised the second half of that.
- **Update (resolved):** follow-up work (see "Recommended follow-up" item
  1) ruled out mesh resolution, both advection schemes, and the
  semi-implicit pressure-velocity coupling as numerical causes of the
  "flat profile" -- the solution converges cleanly under spatial *and*
  temporal refinement. That refinement work turned out to be looking for
  a numerical bug that wasn't there: directly deriving and measuring the
  disturbance's growth rate (item 1, Update 5) shows the solver *does*
  grow the disturbance fast, exactly where the model's own inviscid
  Kelvin-Helmholtz criterion (F > 1, not the different, smaller F0 = 0.36
  friction-enhancement threshold used in the original root-cause finding
  above) says it should -- and correctly stops once the film drains to a
  thickness where that criterion is no longer met. The "flat downstream
  profile" was the converged, physically-consistent behaviour of this
  closure set for that specific case's liquid rate, not an
  under-prediction. Whether the *same* eL-dependent F0-vs-F=1 mechanism
  explains the other 23 originally-failing annular/slug cases is untested
  and a natural next step (item 1, Update 5, final paragraph) -- so the
  classifier's overall accuracy numbers above still stand as measured; what
  changed is the explanation for the representative case, not the score.

## Limitations of this validation pass itself

- Pipe lengths (60D fast sweep, 250D long-pipe subset) are still shorter
  than Shoham's actual test section and far shorter than the source
  paper's own 30 m / D=80mm (L/D ~ 375) validation case; some genuine
  under-development at the fast-sweep length is likely on top of the
  root-cause finding above.
- The regime classifier here uses gas-continuous-fraction (`eg`)
  fluctuation, which is appropriate for the paper's original near-
  horizontal, moderate-holdup validation cases (Section 5) but is
  numerically insensitive for the very thin liquid films / eg close to 1
  that dominate the annular-labeled points in this dataset (an eg bounded
  in [0.997, 1] has little room to show a large absolute fluctuation even
  when the underlying liquid film height varies several-fold). A liquid-
  height-based or relative fluctuation criterion would likely be a fairer
  comparison for those cases and is a good next step alongside the root
  cause above.
- Only flow-*pattern labels* are validated here, not quantitative fields.
  **Update:** see the new section below -- a separate quantitative
  validation pass now compares predicted liquid holdup and frictional
  pressure gradient directly against three independent experimental
  datasets.

## Quantitative validation: liquid holdup and pressure gradient

A second, independent validation pass (`validation/validate_quantitative.cpp`)
compares the solver's predicted liquid holdup and frictional pressure
gradient directly against measured values, rather than flow-pattern
labels, from three real datasets obtained after this document's original
"identified but not accessible" note above (the Andritsos & Hanratty and
Kowalski campaigns remain inaccessible in machine-readable form; see the
paper's Section 7 for that search):

1. **Kokal (1987)**, PhD thesis, University of Calgary -- 1809 runs,
   3 pipe diameters x 7 inclinations, digitized from the thesis via OCR
   (72 out-of-range rows excluded, never "corrected"). 168-case Vsg-spaced
   subsample used.
2. **Newton (1997)**, PhD thesis, UNSW, Appendix A -- 55 horizontal-only
   runs in 50/80mm pipes, all used.
3. **Abdul-Majeed (2022)**, Mendeley Data `10.17632/wyfdm5ysh6.1` -- a
   compiled slug-holdup database from 21 independent studies (Kokal's own
   contribution excluded to avoid double-counting). 126-case subsample.

All 349 cases (`validation/quant_unified_cases.csv` ->
`validation/quant_results.csv`) completed with **zero numerical
failures**. Holdup MAE: Kokal 0.090 (R^2=0.81), Newton 0.113, Abdul-Majeed
compilation 0.303. Error grows monotonically with experimental holdup
magnitude (0.06 at low holdup to 0.29 at high holdup, pooling all three
datasets) -- consistent with, and an independent quantitative confirmation
of, this same document's root-cause finding above: the solver is accurate
in the stratified-like, low-holdup regime and specifically degrades as
the true state moves toward the high-holdup, slug-dominated regime the
flow-pattern classifier was already shown never to identify correctly.
Pressure-gradient agreement is weaker (as expected, since it depends on
holdup through the friction closures): Newton's horizontal cases are
within a factor of 3 in 96% of cases; Kokal's inclined cases, which mix
stratified-like and slug-dominated points at every angle, in 63%. Full
methodology, tables, and figures are in the paper (Section 7,
"Quantitative Validation Against Independent Experimental Campaigns").

## Recommended follow-up

1. Investigate the explicit momentum scheme's linear growth rate for a
   seeded interfacial disturbance against the Kelvin-Helmholtz prediction
   the Andreussi & Persen closure is itself based on, for the specific
   case identified above (D=51mm, Vsl=0.0025, Vsg=25, 1 deg) -- e.g. by
   tracking disturbance amplitude vs. z on a short, well-resolved pipe
   section and comparing the growth rate to the inviscid KH rate implied
   by F. **Update:** we tried the mesh-resolution angle first (see below)
   and it did *not* resolve this, narrowing the likely cause toward the
   scheme's disturbance-growth treatment itself rather than resolution.
   **Update 2:** we then tried MUSCL/TVD flux-limited advection in the
   field continuity equations (`SolverOptions::advectionLimiter`, see
   README.md's "Higher-order / flux-limited advection") on the same stuck
   case. This *did* reveal substantially more holdup-field structure than
   plain upwind (peak gradient up ~50x), but tracking it over time shows
   that structure is a sharper-resolved transient (the solver draining
   from its initial condition toward the true thin-film inlet state) that
   saturates within ~2x rather than continuing to amplify -- not genuine
   growing wave behaviour. That rules out general numerical diffusion in
   the *mass-transport* equations as well, narrowing the cause specifically
   to the **layer momentum equations'** own advection treatment
   (`updateLayerMomentum()`), which was deliberately left untouched that
   pass (non-conservative `u*du/dz` form needs a different TVD derivation
   than the flux-form limiting used for continuity).
   **Update 3:** we then extended the same `advectionLimiter` option to
   the layer momentum equations too, recasting `u*du/dz` into the
   equivalent conservative flux form `d(u^2/2)/dz` (see
   `buildAdvectiveFlux()` in `FourFieldSolver.cpp` and README.md) so the
   same Sweby reconstruction applies. On the stuck case this pushes the
   transient's peak gradient roughly another 5x higher again, but the
   qualitative picture is unchanged: it still saturates (within ~2x)
   rather than continuing to amplify once the transient passes -- no
   sustained wave growth. With both the continuity and layer-momentum
   advection schemes now ruled out as the cause, the two remaining
   untried candidates are the semi-implicit pressure-velocity coupling's
   own damping effect on disturbance growth, and directly measuring the
   seeded disturbance's growth rate against the inviscid KH prediction (the
   original ask in this item) -- which would also clarify whether the
   seeded disturbance's amplitude/frequency is simply a poor match to this
   condition's most-unstable wavelength, independent of the numerical
   scheme entirely.
   **Update 4:** we then investigated the pressure-velocity coupling
   directly, two ways, at the (numerically stable) baseline advection
   scheme. First, a **time-step refinement test**: if the segregated,
   once-per-step (not sub-iterated) pressure correction, or the lagged
   (start-of-step) interfacial-shear coefficient `computeClosures()`
   feeds into that same step's momentum update, were numerically damping
   growth via a step-size-dependent lag, shrinking dt well below the
   CFL-limited default should recover more of the true growth rate.
   Instead, the peak downstream gradient over the run **converges
   cleanly as dt shrinks** (courantTarget 0.5 -> 0.1 -> 0.02, a 25x finer
   step: peak gradient 0.01440 -> 0.01427 -> 0.01426 /m, agreeing to
   <0.1%) -- the time-marching scheme is behaving as a *consistent*
   approximation to whatever the underlying (spatially discretized)
   system actually does, not introducing a growth-suppressing artifact
   that fades as dt -> 0. Second, a **direct `pressureRelaxation` /
   `maxPressureChangeFraction` sweep** at fixed baseline dt: removing the
   under-relaxation entirely (`pressureRelaxation`: 0.2 -> 1.0, i.e. a
   full, unrelaxed correction every step) raised the peak gradient only
   ~19% (0.0153 -> 0.0183 /m) -- a small, bounded, monotonic effect, not
   the difference between "stuck" and genuinely growing; the per-step
   `maxPressureChangeFraction` cap had *no* measurable effect at any
   tested value (0.05, 0.2, 1.0 all gave identical results), meaning it
   was never actually binding at this condition. **This rules out the
   semi-implicit pressure-velocity coupling as the cause too.**
   With mesh resolution, both advection schemes, and the pressure-velocity
   coupling all now ruled out by direct numerical experiment -- and the
   solution behaving as a *converged*, scheme-independent result under
   both spatial and temporal refinement -- the balance of evidence has
   shifted: this increasingly looks like a property of the *converged*
   solution these closures actually produce at this condition (the
   Andreussi & Persen `F` correlation is an empirical friction-factor
   enhancement, fit to steady wavy-stratified data, not a first-principles
   linear-stability growth rate -- `F > F0` need not imply this closure
   combination's converged solution amplifies a small perturbation
   exponentially) or a poor match between the seeded disturbance's
   amplitude/frequency and whatever this system's actual most-unstable
   mode is, rather than a fixable numerical-scheme problem. Directly
   measuring the disturbance's growth rate against the inviscid KH
   prediction (the original ask in this item) is now the most direct way
   to distinguish those two remaining explanations, and is the specific
   next step.
   **Update 5 (resolves this item):** we derived a growth rate to compare
   against, then measured it directly.
   - *Derivation.* Linearizing this model's own layer continuity +
     momentum equations (`updateContinuity()`/`updateLayerMomentum()`),
     in the inviscid limit (no wall/interfacial friction, no entrainment
     -- the same limit the `F` closure's own threshold is based on),
     around a uniform stratified base state gives a 2-layer dispersion
     relation for perturbations `exp(i*k*z - i*omega*t)`:
     `(a+b)*omega^2 - 2k(a*u1+b*u2)*omega + [k^2(a*u1^2+b*u2^2) - C] = 0`,
     with `a = rho_l/A1`, `b = rho_g/A2`, `C = (rho_l-rho_g)*g*cos(theta)*k^2/Si`.
     The instability threshold this implies, `(Ug-Ul)^2 > (rho_l-rho_g)*
     g*cos(theta)*(A2/rho_g + A1/rho_l)/Si`, is k-independent (as it must
     be, matching `F`'s own k-independence) and reduces *exactly* to the
     code's own `kelvinHelmholtzParameterF() > 1` criterion in the
     heavy-liquid limit (`rho_l/A1 >> rho_g/A2`) -- confirmed numerically
     at the stuck case's conditions, where that ratio is ~16,000, i.e. an
     excellent approximation there. This is a useful cross-check that the
     derivation is self-consistent with the code it's meant to describe,
     not an independent, unrelated formula. Where unstable, the growth
     rate is `sigma = Im(omega) = sqrt(4*(a*b*k^2*(u1-u2)^2 - (a+b)*C)) /
     (2*(a+b))`.
   - **Key finding: `sigma(k)` is unbounded, growing linearly with `k`**
     (shorter wavelengths always grow faster, with no most-unstable finite
     wavelength) -- this model's governing equations, as coded, have no
     surface-tension or other short-wave-regularizing term (the paper's
     own Appendix A drops the slip-flux term as small, and no interfacial-
     pressure-difference term is present either), so the continuum problem
     is the classic **short-wave ill-posed two-fluid model** documented in
     the literature (Stewart & Wendroff 1984; Ramshaw & Trapp 1978) for
     exactly this reason: without such a term, growth rate diverges as
     wavelength shrinks toward zero.
   - **`F` at the case's own conditions is above the *true* inviscid
     threshold, not just F0.** Querying `kelvinHelmholtzParameterF()`
     directly at the inlet condition (eL = 0.05) gives **F = 1.14**,
     genuinely past F = 1 -- this is a stronger statement than the
     original root-cause finding's "F = 0.66 > F0 = 0.36", which used a
     downstream, already-drained eL and the wrong (empirical, friction-
     enhancement) threshold for judging exponential growth.
   - **Direct numerical confirmation the solver does grow disturbances
     fast when the model says it should.** Seeding a small perturbation
     (as an initial condition, `bc.seedDisturbance=false` so no inlet
     forcing confounds the measurement) at the same base state, on a
     well-resolved mesh (200 cells/wavelength, `lambda=20D`), the default
     (first-order upwind) solver amplifies it explosively: even starting
     from *pure floating-point roundoff* (~1e-13, no deliberate seed at
     all) it reaches macroscopic amplitude within about 1-1.5 s, at both
     that fine resolution and the original 1-cell/diameter resolution --
     confirming growth is genuinely present in the solver's dynamics at
     this condition, not suppressed by either scheme or mesh. (The
     measured early-time rate, ~1.2-1.4 /s and accelerating, ran faster
     than the single-wavelength `sigma(20D)=0.62/s` prediction and kept
     accelerating over time, consistent with shorter, faster-growing
     content increasingly dominating, per the ill-posedness above -- so
     this confirms growth *happens* and *fast*, without claiming a clean
     single-mode quantitative match, which the ill-posedness makes
     impossible to isolate cleanly by construction.)
   - **The resolution: F is not constant along the pipe -- it falls as
     the film drains, and crosses back below the true threshold.**
     Computing F(eL) at this case's fixed Vsl/Vsg/theta: F = 1.14 at
     eL = 0.05 (the arbitrary initial/inlet condition), crossing **F = 1
     at eL ~ 0.027**, and continuing down to **F = 0.63 at eL = 0.0025**
     (matching the original root-cause finding's downstream F = 0.66
     almost exactly -- confirming this is the same regime). Because
     Vsl = 0.0025 m/s is small, the film's true equilibrium thickness is
     far below the arbitrary eL = 0.05 starting point, and as it drains
     toward that equilibrium -- the same "draining transient" observed
     throughout the mesh/advection/pressure-coupling follow-up work above
     -- F drops below the true instability threshold before the film
     reaches the thickness the original finding queried. **The solver
     was never under-predicting growth: it grows the disturbance while
     F > 1 near the inlet, then correctly stops once the draining film's
     own F drops below 1**, which happens well before it reaches the
     point originally used to diagnose a "stuck" model. F0 = 0.36 (an
     empirical onset of *enhanced interfacial friction/roughness*, not of
     exponential growth) and F = 1 (the *actual* inviscid growth
     threshold, recoverable from this model's own equations) are two
     different thresholds that happened to get conflated in the original
     framing.
   - **Scope of this resolution.** This explains the specific
     representative case traced through this file in detail. Whether the
     same eL-dependent F-crossing mechanism explains the other 23
     originally-failing annular/slug cases in the Shoham sweep (some may
     have inlet/equilibrium conditions where F stays above 1 throughout,
     where genuine sustained growth -- and a real classifier-vs-solver gap
     -- might still be expected) is untested and the natural next step,
     rather than assuming this generalizes without checking.
   - **Update 6 (a separate, related gap the Update 5 derivation
     surfaced, not part of resolving the stuck case itself):** the
     linearization in Update 5 also showed `sigma(k)` is *unbounded*,
     increasing without limit as wavelength shrinks -- this model's layer
     equations carry no interfacial-curvature/surface-tension term, so the
     continuum problem is short-wave ill-posed (Stewart & Wendroff 1984;
     Ramshaw & Trapp 1978), independent of the stuck case's own F<1
     resolution above. `SolverOptions::enableSurfaceTension` adds an
     opt-in (default off) regularization for this. The first, physically
     literal attempt -- an interfacial pressure jump entering the gas
     momentum equation as a third derivative of the liquid height -- was
     implemented and then abandoned: direct von Neumann analysis (and a
     confirming blow-up in testing that a 240x smaller dt only delayed)
     showed every one-sided discretization tried is unconditionally
     unstable under this solver's explicit time-stepping, a genuine
     pitfall of dispersive terms that upwind bias does not fix the way it
     does for advection. What's implemented instead is a biharmonic
     ("hyperdiffusion") proxy added to the liquid-holdup continuity
     equation, with real (non-positive) eigenvalues everywhere and a
     standard, safe explicit stability bound. It is confirmed stable at
     every resolution tested (including the ones that broke the rejected
     third-derivative attempt) and confirmed not to change default
     behaviour; its unit-tested operator and the linear theory behind it
     are sound, but a controlled seeded-wavelength test at the most
     grid-scale resolutions (3-5 cells/wavelength) did not yet show the
     clean, monotonic short-wave suppression the theory predicts at
     reasonable coefficient values. See README.md, "Optional short-wave
     regularization", for the full account and current status -- treat it
     today as a safe, theoretically-motivated, off-by-default option, not
     yet a validated fix for anything.
2. Try a liquid-height-based (or relative) fluctuation criterion in the
   classifier for thin-film conditions, per the limitation noted above.
3. **Tried, then resolved (see item 1, Updates 2-5):** flux-limited,
   less-diffusive advection was added for both the field continuity
   equations and the layer momentum equations' own advection term, and
   the semi-implicit pressure-velocity coupling was tested directly via
   time-step refinement and a relaxation sweep. None of these changed the
   stuck case's qualitative behaviour, because none of them were the
   cause: directly deriving and measuring the growth rate (item 1, Update
   5) shows the solver already grows the disturbance correctly while
   F > 1, and correctly stops once the draining film's F drops below 1 --
   a real, physically-driven transition in the base state, not a
   numerical deficiency in any of the schemes tested.
4. For a quantitative (not just regime-label) inclined validation,
   transcribe a small set of published holdup/pressure-drop points from
   Barnea, Shoham, Taitel & Dukler's inclined-pipe papers (no open dataset
   found) and extend `validate_shoham.cpp`'s comparison logic to numeric
   fields instead of (or alongside) regime labels.
5. **Tried:** two mesh-adaptivity strategies were added and tested
   directly against the D=51mm/Vsl=0.0025/Vsg=25/1deg case above -- an
   h-refinement AMR (`SolverOptions::amr`) and a moving/r-adaptive mesh
   (`SolverOptions::movingMesh`), both driven by the same Kelvin-Helmholtz
   `F` indicator this section discusses. Neither produced genuine wave
   growth on that case, even the moving mesh at 1-diameter-mesh-equivalent
   fine resolution: the liquid film stayed essentially flat downstream
   regardless of how finely it was resolved. On a *different* case where a
   front does form (the horizontal slug-formation demo), the moving mesh
   demonstrably sharpens the captured gradient by ~8.7x over the fixed
   grid, showing the remap/relocation machinery itself works correctly --
   so the stuck case's flat profile is not a meshing bug. This is a
   meaningful (negative) result: it rules out mesh resolution as the
   explanation for that specific case and points toward the momentum
   scheme's growth-rate treatment (item 1) as the more likely cause. See
   README.md, "Adaptive mesh refinement" and "Moving-mesh (r-adaptive)
   tracking", for the full writeup and how to reproduce both tests.
6. **Done:** implemented `SolverOptions::enableTurbulentViscosity`, a
   physics-based nonlinear regularization adapted from Lopez-de-Bertodano
   & Clausse (2026), "Nonlinear Stability in the Two-Fluid Model of
   Two-Phase Flow" (arXiv:2509.04679) -- see README.md, "Turbulent-
   viscosity nonlinear regularization" for the full derivation,
   translation to this solver's equations, and honest test results.
   Important scoping note carried over from item 1 above: this option
   targets NONLINEAR wave-growth saturation (does a large-amplitude wave
   stay bounded or blow up), not the LINEAR question this section's own
   stuck case already resolved (that case's F<1 base state is genuinely
   stable, so no regularization of this kind can or should "unstick" it).
   Confirmed stable (zero NaN/Inf) across mixing lengths 0-20% of pipe
   diameter on the horizontal slug-formation demo, with a physically
   sensible damping trend.
   **Update:** ran the sustained-F>1 saturation test. First confirmed,
   using `kelvinHelmholtzParameterF()` directly (the same closure the
   solver uses) swept over eL at this case's fixed Vsl/Vsg, that F here
   crosses 1 only once the local holdup climbs to eL ~ 0.75-0.8 -- the
   opposite profile from the stuck case above (which *drains* toward
   stability): here, growth toward a near-blocking holdup is genuinely
   self-accelerating once past that threshold, exactly the mechanism
   nonlinear regularization is meant to address. Ran baseline (no
   regularization), surface-tension-only, and surface-tension +
   turbulent-viscosity (mixing length 10% and 20% of D) out to 45s
   simulated time (order 30,000-35,000 steps each) on this case.
   **Result: all four configurations, including the completely
   unregularized baseline, remain numerically stable throughout**, with
   peak holdup repeatedly climbing to near-full blockage (eL approx.
   0.95-0.999, i.e. a slug-like spike) and then relaxing back to
   0.6-0.7 in a bounded, recurring cycle -- never a monotonic runaway.
   This means the case does *not* cleanly isolate turbulent viscosity's
   contribution: this codebase's existing robustness safeguards
   (velocity clamping, momentum-fraction floors -- see "Robustness
   safeguards" in README.md) already prevent literal numerical blow-up
   on this case regardless of the new physics term, unlike the reference
   paper's bare toy model, which has no such independent safety net and
   genuinely does blow up without it. The turbulent-viscosity
   configurations' peak-holdup trajectories differ modestly from the
   baseline's (generally similar magnitude, not dramatically suppressed)
   at the mixing lengths tested, consistent with a real but modest
   contribution layered on top of an already-bounded system, rather than
   being the sole thing standing between stability and blow-up the way
   it is in the reference paper. A cleaner isolation of this term's own
   contribution would need a case (or a temporarily disabled safety
   clamp) where the unregularized baseline *does* blow up on its own --
   not attempted here, flagged as the honest next step rather than
   claimed.
7. **Done:** implemented `SolverOptions::enableImplicitFriction`, an IMEX
   (implicit-explicit) time-integration split treating wall and
   interfacial friction implicitly (a local 2x2 linear solve per face,
   linearising each nonlinear friction force by freezing its |velocity|
   factor at the old time level) while advection stays explicit -- see
   README.md, "Implicit friction / IMEX time integration" for the full
   derivation and test results. Direct dt-sweep test on a high-drag case
   confirms the intended effect cleanly: at dt just above the case's own
   `stableTimeStep()` (which has no friction-stiffness cap today),
   explicit friction blows up while implicit friction stays stable and
   converges to the same answer smaller, unambiguously-stable dt values
   reach. On typical already-stable cases (e.g. the horizontal
   slug-formation demo) the effect is modest, not dramatic -- comparable
   step count and wall-clock time, physics within 1-3% -- so this is a
   robustness/stability-margin improvement for the stiff-friction regime
   specifically, not a general speedup. A genuine, real bug was caught
   and fixed before relying on any of this: an early draft's linear
   coefficient divided the *signed* friction force by |velocity| instead
   of the force's own magnitude by |velocity|, which would have silently
   flipped the sign of friction (into anti-friction) for negative
   velocities -- found by reading `wallShearStress()`'s and
   `interfacialShearStress()`'s actual implementations rather than
   assuming their sign convention.
8. **Done:** parallelized the outer per-case loop of both validation
   drivers (`validate_shoham`, `validate_quantitative`) with OpenMP
   (`#pragma omp parallel for`, dynamic scheduling) -- see README.md,
   "Parallelization (OpenMP)" for the full rationale, including why the
   *inside* of a single `FourFieldSolver` run was judged not worth
   parallelizing at this project's typical resolution (N=40-300 cells),
   in contrast to the hundreds of fully-independent solver instances each
   validation sweep runs. Confirmed by direct test: full 349-case and
   182-case sweeps produce byte-identical output CSVs (including the
   confusion matrix and accuracy tally) whether run with
   `OMP_NUM_THREADS=1` or `4`, since the outer loop has no shared mutable
   state between cases and results are written to a preallocated
   per-case slot then serialized to disk in a second, purely sequential
   pass. Measured wall-clock speedup on 4 physical cores: 4.0x on
   `validate_quantitative` (294.8s -> 73.8s) and 3.9x on `validate_shoham`
   (55.7s -> 14.2s) -- both close to the ideal 4x, as expected for this
   embarrassingly-parallel a workload.
9. **Done:** implemented `SolverOptions::enableWellBalancedGravity`, a
   discrete-hydrostatic-reconstruction fix for the interface-slope gravity
   term's discretization in terrain-following (spatially varying `theta`)
   runs -- see README.md, "Well-balanced gravity discretization" for the
   full derivation and test results. A standalone equilibrium test (true
   per-cell static-equilibrium profile for a 4-segment, uneven-angle
   terrain pipe, zero net flow, interior-face velocity tracked from rest)
   confirmed the default scheme's spurious drift shrinks at close to first
   order under mesh refinement (1.36e-5 m/s at N=40 down to 1.04e-6 m/s at
   N=320) -- real, but ordinary truncation error, not a structural defect.
   The fix reduces this by roughly five orders of magnitude (6.2e-11 to
   8.3e-11 m/s across the same refinement sweep) and, notably, the residual
   stops shrinking with resolution at that point -- consistent with having
   reached the floating-point noise floor rather than a resolution-limited
   error, the expected signature of *exact* (not just asymptotic)
   well-balancing. Regression-confirmed bit-for-bit identical default (off)
   behaviour against the full demo suite. On the existing terrain demo
   case (a genuinely flowing case that develops a KH-unstable slug front,
   not an equilibrium case), enabling the option left total liquid mass
   within 0.5% and the bulk of the holdup profile closely matched, but
   shifted the exact slug-front position -- diagnosed, via a side-by-side
   profile/mass comparison (not just asserted), as the same
   chaotic-sensitivity signature already documented in item 7's
   floating-point-reordering finding, not a new bug.
10. **Done:** implemented `stepImplicitPressureVelocity[Adaptive]()`, a
    Jacobian-Free Newton-Krylov (JFNK) fully-implicit alternative to
    `step()`'s segregated coupling, directly addressing the limitation
    flagged in "Recommended follow-up" item 4 above and in README.md's
    "Pressure-velocity coupling investigation" -- see README.md,
    "Jacobian-Free Newton-Krylov fully-implicit coupling" for the full
    formulation, the debugging narrative (several plausible causes --
    a GMRES bug, friction-correlation branch discontinuities, saddle-point
    ill-conditioning needing a preconditioner -- investigated by direct
    test and ruled out, before a dt-sweep diagnostic correctly identified
    ordinary Newton-basin shrinkage, not a defect, as the actual cause),
    and the measured results. The decisive test directly replicates the
    exact stiff case and dt (D=0.05 m, L=5 m, N=50, Vsl=0.05, Vsg=20,
    dt=0.005 s, 100 steps) that item 4's own earlier investigation showed
    the explicit scheme blowing up on: the explicit scheme again blows up,
    at step 64/100, while the new adaptive JFNK wrapper completes all 100
    steps stably (`|u2-u1|_max` = 21.72 m/s), at an honestly-reported cost
    of ~56 s wall-clock (vs. the explicit scheme's sub-5 ms before
    diverging) and an average of 61.6 sub-steps / 33.1 Newton iterations /
    1958.6 GMRES matrix-vector products per requested macro-step. On a
    milder, non-stiff case the two schemes' trajectories agree closely at
    matched dt, confirming JFNK converges to the same physics, not merely
    to *a* stable answer. Regression-confirmed bit-for-bit identical
    `step()` behaviour (the new code is isolated in its own translation
    unit and only runs when explicitly invoked).
11. **Done:** implemented `stepSegregatedAccelerated()`, testing whether a
    much cheaper fixed-point accelerator (Richardson/Picard iteration,
    optionally with Anderson mixing) could recover JFNK's stability
    benefit on the same decisive stiff case without a Jacobian-vector
    product or GMRES -- see README.md, "Anderson-accelerated Picard
    solve" for the full account, including two real bugs caught by direct
    testing (a sign error in the Anderson update, found because it made
    Anderson perform *worse* than plain Picard at every depth tried; a
    mis-scaled regularizer in its small least-squares solve) before
    trusting any result built on the code. With both fixed: undamped
    Richardson diverges outright (confirming the diagonal scaling that
    makes JFNK's Jacobian-vector products well-posed does NOT make the
    identity a usable stand-in for the true Jacobian); heavy damping
    (`beta ~ 0.003`) stabilizes it; and Anderson mixing's effect, once
    correctly implemented, is genuine but inconsistent -- helps markedly
    on a mild case at a shallow window, hurts at every depth tested on
    the stiff case that matters, reported as a real negative result
    rather than tuned away. The actual headline result came from the
    PLAIN damped baseline (`andersonDepth = 0`), not from Anderson mixing
    at all: on the identical stiff case/dt/step-count JFNK was validated
    against, it converges all 100 steps to the SAME `|u2-u1|_max = 21.72`
    m/s JFNK-adaptive found, in 1.07s wall-clock -- roughly 50x cheaper
    than JFNK-adaptive's ~56s for an identical outcome. Explicitly caveated
    in README.md as a single-case result (the damping factor and
    tolerance were found by direct sweep on this one case, not derived
    generally) rather than a claim that JFNK's more robust, automatically
    adaptive machinery is unnecessary in general. Regression-confirmed
    bit-for-bit identical `step()` behaviour throughout.
12. **Done:** implemented `SolverOptions::enableETDFriction` (exact
    matrix-exponential integration) and `enableIMEXRKFriction` (a 2-stage,
    2nd-order, L-stable SDIRK, an instance of the Pareschi & Russo 2005
    IMEX-RK framework) as two further alternatives to
    `enableImplicitFriction`'s single-stage backward Euler, all three
    sharing the identical friction linearization -- see README.md,
    "Exponential time differencing and IMEX Runge-Kutta friction" for the
    full account. Caught a genuine bug in the first `enableIMEXRKFriction`
    attempt (a plausible-looking but non-L-stable Butcher tableau) not by
    deriving its stability function up front, but because the stiff dt-
    sweep test below showed it blowing up where plain backward Euler
    stays stable -- the wrong direction for a strictly more accurate
    integrator of the same linear system, which prompted the derivation
    that found the actual tableau error; the corrected scheme was then
    verified L-stable both analytically and by a direct numerical sweep
    before being trusted. Once both were correct, found and reported
    honestly a genuine negative result rather than a win: on a mild case
    all three friction integrators show comparable error against a fine-
    dt reference (the dominant error there is the shared linearization,
    not the choice of integrator); on the established stiff benchmark
    case, step-by-step tracking of `|u2-u1|_max` (not just pass/fail)
    shows backward Euler locking onto a stable fixed point from step 0,
    while plain explicit, ETD, and IMEX-RK all instead show the same
    quantity growing slowly over dozens of steps before eventually
    diverging -- diagnosed as backward Euler's own lower formal accuracy
    acting as unintentional extra numerical dissipation that happens to
    suppress a genuine, slowly-growing instability in this strongly-
    sheared, low-holdup regime (the same regime item 1's own
    Kelvin-Helmholtz case study examines at length), which the more
    accurate integrators correctly do not suppress. `enableImplicitFriction`
    remains the recommended choice for this solver's stiff operating
    envelope; ETD and IMEX-RK are documented as a well-tested negative
    result and as independently-verified, reusable closed-form building
    blocks, not as replacements. Regression-confirmed bit-for-bit
    identical `step()` behaviour throughout.
13. **Done:** implemented `SolverOptions::enableDataDrivenClosureCorrection`,
    rescaling only the `AndreussiPersen1987` interfacial friction
    correlation's enhancement above the baseline gas-wall friction factor
    -- see README.md, "Data-driven closure correction" for the full
    account. Fit by a genuine train/test line search: trained on Kokal
    1987 + Newton 1997 (223 cases), evaluated on the fully independent,
    never-fit-against 126-case Mendeley set. Result is a clean negative
    finding, reported as such rather than tuned until it looked like a
    win: the training-set response to the scale parameter is nearly flat
    and non-monotonic across a 10x range (0.0 to 3.0), and the value that
    scores best on the training data (2.0, MAE 0.0944 vs. the unmodified
    closure's 0.0960) makes the held-out test set very slightly *worse*
    (0.3047 vs. 0.3032 unmodified), while a value that scores *worse* on
    training (0.0, removing the enhancement entirely) gives the lowest
    test error of anything tried (0.2887) -- opposite rankings on the two
    sets, the classic signature of overfitting a single scalar to noise,
    not a real physical relationship. Conclusion: the stratified-vs-slug
    accuracy gap the "Quantitative validation" section above already
    identified is not primarily attributable to this correlation's
    enhancement *magnitude*
    being miscalibrated by a single multiplicative factor; closing it
    would need something more structural (a different functional form,
    a correction to the F0=0.36 threshold itself, or a fix to an
    entirely different closure). `closureCorrectionScale`'s default
    stays at 1.0 (the unmodified closure) because no tested value
    actually generalized. Regression-confirmed bit-for-bit identical
    `step()` behaviour throughout.
14. **Done:** tried the F0=0.36 threshold correction item 13 flagged as
    a possible next step -- `SolverOptions::closureF0Shift`, shifting
    the AndreussiPersen1987 correlation's own Kelvin-Helmholtz onset
    threshold (`F > F0 + closureF0Shift` instead of `F > F0`) rather
    than rescaling the enhancement's size -- see README.md, "A second
    attempt: shifting the F0 threshold itself" for the full account.
    Same train (Kokal+Newton, 223 cases) / held-out test (Mendeley, 126
    cases) methodology as item 13; same kind of result. The train-set
    response to the shift is just as flat as the magnitude knob's (MAE
    0.0954-0.0970 across a `-0.3` to `+1.0` sweep, barely distinguishable
    from the unmodified closure's 0.0960), and the nominal best value
    (`-0.05`) makes the held-out test set very slightly worse (0.3036 vs.
    0.3032 unmodified) -- the same overfitting signature item 13 already
    found, now confirmed on an independent parameter of the same closure.
    Conclusion, now doubly confirmed rather than a single inconclusive
    trial: the stratified-vs-slug accuracy gap is not hiding in how the
    AndreussiPersen1987 correlation is calibrated, in either of its two
    natural free parameters; closing it would need a different closure
    entirely, a structural change to the model, or scrutiny of the
    validation driver's own holdup sampling -- not a scalar tweak to
    this one correlation. `closureF0Shift`'s default stays at 0.0 for
    the same reason. Regression-confirmed bit-for-bit identical
    `step()` behaviour throughout.
