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

**Root-cause finding.** A representative failing case (D = 51 mm, Vsl =
0.0025 m/s, Vsg = 25 m/s, 1 deg incline, labeled annular) was inspected
directly: even at 1-diameter mesh spacing and 10+ pipe residence times,
the liquid film height along the back half of the pipe is essentially flat
(1-2% variation -- numerical noise, not a growing wave). Yet querying the
solver's own interfacial-friction closure at the same local conditions
gives F = 0.66 against the Andreussi & Persen (1987) threshold F0 = 0.36
(see `Closures.hpp`), i.e. **the closure itself reports this point as
being in the enhanced-shear, unstable regime** (an 11.9x enhancement of
the interfacial friction factor over its smooth-wall value) -- but the
explicit time-marching solution is not amplifying the seeded disturbance
into a growing wave. That combination (unstable-per-closure, stable-in-
practice) points to the disturbance not growing fast enough under this
codebase's explicit/semi-implicit momentum treatment, rather than to a
wrong closure choice. This is a genuine, reproducible limitation of the
current implementation, not a data or methodology artifact, and is the
most concrete lead for follow-up work (see "Recommended follow-up" below).

## Interpretation

- The **flow-regime-independent stratified prediction is strongly
  validated across the entire inclination range this codebase adds**
  (100% on horizontal, 85% on inclined, using the *same* closures at every
  angle, no regime-specific tuning) -- directly supporting the paper's
  central claim and the specific horizontal/inclined comparison asked for.
- **Wave growth into annular/slug/bubbly regimes is under-predicted** in
  the current implementation once conditions call for it, at every
  inclination (not specifically worse at steep angles) -- this is a real,
  now well-localised gap, not a vague "needs more validation" caveat.
- The picture is therefore: *the model correctly stays stratified when
  the flow should be stratified, and correctly resolves stratified
  behaviour identically well from horizontal through steep incline, but
  currently under-predicts the growth of disturbances into the other
  regimes.* That's a meaningfully different (and more actionable)
  conclusion than either "it works" or "it doesn't".
- **Update:** follow-up work (see "Recommended follow-up" item 1) has
  since ruled out mesh resolution, both advection schemes, and the
  semi-implicit pressure-velocity coupling as the cause of the
  under-predicted growth, each via a direct numerical experiment rather
  than by inspection. The solution converges cleanly under spatial *and*
  temporal refinement, which points away from "the numerical scheme is
  suppressing real growth" and toward this being what these closures'
  converged solution actually does at this condition -- i.e. a modelling/
  closure question now, more than a numerical one.

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
- Only flow-*pattern labels* are validated here, not quantitative fields
  (holdup, pressure drop, slug frequency/length). The source paper's own
  Section 5 already validates several of those quantitatively (Nydal et
  al. 1992 slug length/frequency and Bendiksen 1984 bubble velocity, both
  horizontal); repeating that quantitative comparison for inclined cases
  would need a quantitative inclined dataset, which -- per the search
  above -- exists in the literature (Barnea/Shoham/Taitel/Dukler) but not
  as an open, machine-readable dataset.

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
2. Try a liquid-height-based (or relative) fluctuation criterion in the
   classifier for thin-film conditions, per the limitation noted above.
3. **Tried (see item 1, Updates 2-4):** flux-limited, less-diffusive
   advection was added for both the field continuity equations and the
   layer momentum equations' own advection term, and the semi-implicit
   pressure-velocity coupling was tested directly via time-step
   refinement and a relaxation sweep. None of these produced genuine wave
   growth on the stuck case, or even meaningfully changed its magnitude
   -- ruling out numerical diffusion in the advection schemes *and* the
   pressure-velocity coupling as the explanation. See item 1's Update 4
   for what that leaves as the likely remaining explanation.
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
