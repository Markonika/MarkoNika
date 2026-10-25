# Flat 2D planar approximation — plan (not implemented)

This is the fallback plan if the thin 3D wedge case turns out to be too
expensive to iterate on quickly (large mesh, slow turbulence-model
convergence, etc.). It is **not built** — this file is a plan to execute
later if needed, kept separate from the actual 3D wedge case in this
directory (`msfr_wedge.flml`, `generate_mesh.py`).

## What "flat 2D planar" means here, precisely

Take the exact same (r, z) profile curve already built for the wedge
(`generate_mesh.py`'s `smooth_T(r)` function and the duct geometry), but
instead of revolving it around the vertical reactor axis, mesh it as an
ordinary flat 2D domain — literally substituting Cartesian (x, y) for
(r, z) and ignoring that it was ever meant to represent a radial cross-
section of a revolved body.

**This is a real physics simplification, not just a cheaper mesh.** The
true governing equations in the wedge case include the 1/r geometric
spreading/contracting of streamlines as fluid moves to larger or smaller
radius (conservation of mass in an annular control volume scales with
`2*pi*r*dr`, not `dr`). A flat 2D domain has no such effect baked in —
it's dimensionally and physically a different flow (effectively an
infinite-depth planar duct, not a torus). Concretely, this means:
  - Near the inner wall (smallest r), the real axisymmetric flow
    *accelerates* as it's squeezed toward the axis; the flat 2D version
    won't show this.
  - The flat 2D version cannot be directly compared, even qualitatively
    in the near-axis region, to the velocity/temperature fields reported
    in Rouch et al. (2014) or Fang et al. (2024), which are true r-z
    axisymmetric (or 3D wedge/full-core) results.
  - It is still useful as: (a) a much cheaper way to validate mesh
    generation, boundary condition wiring, and the anisotropic
    adaptivity workflow in Fluidity before committing to the 3D wedge
    case, and (b) a qualitative (not quantitative) demonstration of
    adaptivity resolving the same recirculation/shear-layer topology,
    since the overall shape (duct in, big curved cavity, duct out) is
    preserved.

## Why this needs a separate mesh (not just "drop the revolve step")

`generate_mesh.py`'s profile curve is built directly in a plane (the
`add_pt` function assembles `(r, z)` points, which then get embedded as
3D coordinates `(z, r, 0)` and revolved). For the flat 2D plan, the
fix is simpler: skip the revolve entirely and mesh the planar surface
(`surf` in the script, before the `geo.revolve(...)` call) directly as a
2D mesh, with physical *lines* (not surfaces) for inlet/outlet/walls,
and no symmetry-plane boundaries needed at all (there's no "thin wedge"
in 2D — the 2D domain has only 4 boundary curve groups: inlet, outlet,
wall-top, wall-bottom, which can be lumped into a single "wall" group
like the 3D case does).

## Steps to build it (left for later; summarized, not executed)

1. In a copy of `generate_mesh.py`, delete everything from the
   `# Revolve about the x-axis...` comment onward. Keep the profile
   construction identical (same `R_MIN`, `R_MAX`, `H_IN`, `H_OUT`,
   `DUCT_HEIGHT`, `DUCT_LENGTH` — same sourced envelope and the same
   placeholder duct dimensions, so results stay comparable across the
   two cases).
2. Tag physical *curves* (1D, not surfaces) for `inlet`, `outlet`, and
   `wall` directly on the loop's line entities (`l_inner_cap`,
   `l_top_wall`, ..., `l_bot_wall`), and tag the plane `surf` itself as
   the single physical 2D `fluid` region.
3. Mesh with `gmsh.model.mesh.generate(2)`, write `.msh` v2.2 as before.
4. Build a `.flml` nearly identical to `msfr_wedge.flml`, but:
   - `<dimension>` → 2
   - Drop the `SymmetryPlanes` boundary condition block entirely (no
     equivalent in 2D — there's nothing to be symmetric about).
   - Velocity initial condition / dirichlet values become 2-component
     (`shape="2"`) instead of 3-component.
   - The inlet parabolic-profile Python function simplifies to use only
     `X[0]` (axial) and `X[1]` (the single in-plane "radial-like"
     coordinate) rather than `math.hypot(X[1], X[2])`.
   - Everything else (material properties, placeholder effective
     viscosity, mesh-adaptivity block borrowed from
     `driven_cavity_2d_vms_anisotropic.flml`) carries over unchanged,
     since none of it is inherently 3D.
5. This 2D case is also the natural one to attach `mesh_adaptivity` /
   `hr_adaptivity` to directly (the project's existing anisotropic
   adaptivity work is itself 2D), without needing to extend that
   machinery's 2D-specific code paths to 3D first.

## When to actually build this

Switch to this plan if the 3D wedge case (`msfr_wedge.flml`) turns out
to need mesh sizes or turbulence-model iteration counts that make
quick design iteration impractical, or if the goal shifts toward
exercising the adaptivity machinery itself (which is presently 2D-only
in this codebase) rather than getting a physically faithful MSFR flow
field.
