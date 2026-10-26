# MSFR Nek5000 reference data (INL Virtual Test Bed)

**Source:** INL Virtual Test Bed, "MSFR Turbulence Modeling Parameter Tuning
Using Nek5000 Data"
<https://virtualtestbed.inl.gov/msr/msfr/griffin_pgh_results.html>
(page could not be fetched directly from this sandbox -- `*.inl.gov` is
blocked by this environment's network policy; the user opened it in their
own browser and provided the page text and figures on 2026-10-26).

**Geometry:** same mesh as the Nek5000 CFD model page
(<https://virtualtestbed.inl.gov/msr/msfr/nek5000_cfd_model.html>), i.e. the
real MSFR core cavity geometry -- same source family as
`msfr_wedge/vtb_source/msfr_v2.geo` already in this repo. The plotted domain
is the **2D axisymmetric (r,z) cross-section** (full "bucket" profile shown
mirrored about the reactor's symmetry/revolution axis), not a 3D wedge.

## What this gives us

Two RANS-vs-Nek5000 comparison cases are documented, at Re=40,000 and
Re=1,000,000 (the latter matching the real MSFR operating regime). Boundary
conditions described on the page, Re=40K case:

- velocity inlet profile
- pressure outlet
- symmetry at the center of the geometry along the (revolution) axis
- no-slip at walls (Re=40K); standard velocity wall functions at walls
  (Re=1M)
- time derivatives added to the momentum equation for relaxation
- Re=1M case additionally uses an exponentially decaying viscosity ramp for
  convergence: `mu_sim = mu_phys * (exp(-3*t) + 1)`

This is structurally very close to what our own Fluidity wedge case already
does (velocity inlet, pressure-only outlet, symmetry BC on the cut faces,
no-slip walls) -- the main difference is this Nek5000 case is a true 2D
axisymmetric solve, while our Fluidity case approximates that with a thin 3D
wedge (Fluidity has no native axisymmetric solver).

**δ (mixing-length cap) tuning result:** δ=0.12 m at Re=40,000, δ=0.10 m at
Re=1,000,000, obtained by minimizing Pronghorn's velocity-field error against
Nek5000 specifically at the **minimum core radius** cross-section (the
r_core=1.054 m waist of the geometry).

## Files in this directory

- `nek5000_velocity_magnitude_re40k.png` / `pronghorn_velocity_magnitude_re40k.png`
  -- velocity magnitude field, Re=40K, Nek5000 vs. Pronghorn.
- `nek5000_pressure_re40k.png` / `pronghorn_pressure_re40k.png`
  -- pressure field, Re=40K, Nek5000 vs. Pronghorn.
- `velocity_at_min_core_radius_re40k.png` -- line plot of velocity
  magnitude vs. radius `r`, sampled at the minimum core radius
  cross-section, Re=40K, both codes. This is the only plot on the page
  with an explicit numeric axis (velocity magnitude 0-1.4, r 0-~1.05 m),
  so it's the one quantitative dataset we could extract.
- `velocity_at_min_core_radius_re40k_digitized.csv` -- **digitized from
  the PNG above** by color-matching the Nek5000 (orange, RGB
  255,127,14) and Pronghorn (blue, RGB 31,119,180) curves pixel-by-pixel
  and calibrating against the plot's gridlines (not extracted from any
  underlying numeric data file -- none was published). Columns:
  `r_m`, `velocity_magnitude_pronghorn`, `velocity_magnitude_nek5000`,
  45 points evenly spaced over r in [0, 1.054] m. Treat this as
  approximate (pixel-digitization error, probably a few percent) --
  good enough for a qualitative/semi-quantitative comparison shape, not
  for precision validation.

## Caveats

- No raw field data (`.vtu`/`.vtk`/exodus) was published alongside this
  page for the Re=40K or Re=1M cases -- only these rendered plots. The
  page separately mentions sample outputs for the *coupled Griffin-Pronghorn*
  model (not the Nek5000-vs-Pronghorn tuning study) are hosted on git-lfs,
  which is a different, much larger multiphysics case, not a drop-in
  CFD-only comparison.
- Units/normalization of the velocity and pressure fields (e.g. whether
  pressure is relative/non-dimensional) are not stated on the page as
  transcribed; treat absolute values with caution and prefer comparing
  *shape* (recirculation pattern, peak location) first.
- Re=40,000 is far below the real MSFR operating Reynolds number (~1e6,
  also shown on the source page) but is a much more tractable target for
  a first quantitative comparison against our own Fluidity case, which has
  so far used an unsourced placeholder viscosity.
