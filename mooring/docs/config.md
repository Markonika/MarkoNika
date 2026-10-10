# Configuration reference (JSON, SI units)

Run: `build/mooring_run config.json [key.path=value ...]` (overrides are JSON values, e.g. `motion.radius_m=0.1`).
Outputs (if `output.write` is true): `<tag>_timeseries.csv` (time, top tension, anchor tension, top x/y/z; raw,
decimated to `dt_out_s`), `<tag>_cycles.csv` (per-cycle maxima from every time step), `<tag>_params.json` (full effective
configuration plus derived values, time step, slack/clipping/contact counters).

| Key | Meaning |
|---|---|
| `line.length_m`, `line.segments`, `line.EA_N`, `line.mass_per_length_kg_m` | unstretched length, N, axial stiffness, mass per length |
| `line.density_kg_m3` | cable density -> submerged weight, Eq. 3.26 (else `line.weight_per_length_N_m`) |
| `line.hydro_diameter_m` (D0), `line.soil_diameter_m` (D1), `line.nominal_area_m2` (A1, default pi/4 D0^2) | Morison / soil diameters |
| `line.Cm`, `line.Cdt`, `line.Cdn` | added mass, tangential and normal drag coefficients |
| `line.tension_curve`, `line.maxwell_branches` | optional rope law (below); `EA_N` may be omitted when a curve is given |
| `line.internal_damping_Ns` | c_int in T = EA eps + c_int d(eps)/dt (0 = off; also the slack regularisation, logged) |
| `environment.hydro`, `.seabed` | enable Morison forces + submerged weight; enable seabed |
| `environment.rho_w_kg_m3`, `.water_surface_z_m`, `.seabed_z_m` | water density, still-water level, flat seabed height |
| `soil.stiffness_Pa_per_m`, `.damping_factor`, `.friction`, `.v_lim_m_s` | Ks, zeta_s, mu_s, v_lim (Eqs. 3.32-3.36) |
| `anchor_m`, `fairlead_rest_m` | anchor (fixed) and top end at rest |
| `initial.shape` | `touchdown` (catenary lying on the bed, then static relaxation) |
| `motion.type` | `none`, `circle_xz` or `harmonic` (sinusoidal along `direction_vec`, amplitude `radius_m`) (top end on a circle about `centre_m`, default = rest position) |
| `motion.plane_angle_deg` | rotates the circle plane about the vertical axis (set `fairlead_rest_m` consistently) |
| `environment.current_m_s` | uniform steady current vector [m/s] (drag only; also applied in the static solve) |
| `initial.perturbation_y_m` | out-of-plane half-sine perturbation of the relaxed line, released from rest (3D only) |
| `motion.radius_m`, `.period_s`, `.direction` (+1 counter-clockwise in x-z, -1 clockwise), `.phase_deg`, `.ramp_cycles`, `.cycles` | excitation |
| `numerics.planar` | true: 2D mode, motion confined to the plane y = y_anchor (forces/velocities in y dropped) |
| `numerics.scheme` (`rk4`/`verlet`), `.cfl`, `.dt_s` (override), `.relax_force_tol`, `.relax_max_steps` | integration and static relaxation |
| `point_elements[]` | list of `{type: "floater"\|"clump"\|"generic", node \| arclength_m (nearest node), mass_kg, diameter_m, Cd, Cm, buoyancy_N (floater) \| submerged_weight_N (clump) \| volume_m3, area_m2 (generic)}` |
| `statistics.first_cycle` | first cycle included in the mean of the cycle maxima |
| `output.directory`, `.tag`, `.dt_out_s`, `.write` | output control |

## Platform runner (`build/mooring_platform config.json [key.path=value ...]`)
One file holds the body, the lines and their fairleads. Outputs: `<tag>_timeseries.csv` (time, 6 DOF, and per line the fairlead force
magnitude and the end-segment tension) and `<tag>_params.json` (full config + derived values: dt_body, shortest line step, sub-step
ratio, coupling scheme, equilibrium status).

| Key | Meaning |
|---|---|
| `body.mass_kg`, `body.inertia_diag_kg_m2` (3) or `body.inertia_kg_m2` (9) | rigid-body mass and inertia about the CG (reference axes) |
| `body.cg_ref_m` | global position of the CG in the reference (free-floating equilibrium) pose |
| `body.A`, `.B`, `.C` (6x6 or 36 numbers) or `body.A_diag`, `.B_diag`, `.C_diag` (6) | added mass, damping, hydrostatic stiffness (constants) |
| `body.Dq` (6), `body.F0` (6) | quadratic drag per DOF; constant generalised load |
| `lines[]` | `{name, anchor_m, fairlead_body_m (relative to the CG, reference axes), line: {...}, point_elements: [...], initial: {shape}}` |
| `line_defaults` | `line` fields shared by all legs; each leg's `line` overrides them |
| `environment`, `soil`, `numerics` (top level) | defaults for every line (same keys as the single-line runner) |
| `initial.xi0`, `.xi_dot0` | pose / velocity before the equilibrium solve |
| `initial.equilibrium`, `.equilibrium_tol`, `.xi_offset` | solve the static equilibrium, then displace by `xi_offset` (e.g. for a decay test) |
| `numerics.dt_body_s`, `.t_end_s` | body step and end time; `numerics.cfl`/`line_dt_s`/`scheme` set the line integration |
| `output.directory`, `.tag`, `.dt_out_s`, `.write` | output control |

Example: `examples/platform/three_leg_example.json` (generic, **not** the Paredes buoy).

## Waves (milestone 9, platform runner)
| Key | Meaning |
|---|---|
| `waves.depth_m`, `waves.surface_z_m` | constant water depth and still-water level (default: `environment.water_surface_z_m`) |
| `waves.components[]` | `{height_m` or `amplitude_m, period_s, phase_rad, direction_deg}` Airy components (sum = irregular kinematics) |
| `waves.stretching` | `wheeler` (default) or `none` |
| `waves.ramp_time_s` | cosine ramp of all amplitudes over the first seconds (0 = off) |
| `body.wave_force.w`, `.delta` | 6 + 6 numbers: `f_i = w_i A sin(omega t + delta_i)` with A and omega of the first component (Eq. 3.60) |
| `initial.xi_dot0` | initial body velocity; applied **after** the equilibrium solve as well |
Elevation of a component is `A sin(omega t - k s + phase)` (s along the direction): at the origin A sin(omega t). The water velocity/acceleration feed the lines and point elements.

### Point-element anisotropy keys (generic type)
`Cd_t`, `area_t_m2`, `Cd_n`, `area_n_m2` (drag along/normal to the line; used when `Cd_t` or `Cd_n` is set) and `added_mass_t_kg`, `added_mass_n_kg` (anisotropic added mass). Absent keys keep the isotropic `Cd`/`area_m2`/`added_mass_kg` behaviour.

### Rope law (nonlinear / viscoelastic axial behaviour)
`line.tension_curve`: `[[strain, tension_N], ...]`, strictly increasing and positive, piecewise linear through (0, 0), last slope extrapolated; replaces `EA_N` (if both are given `EA_N` is only the reference stiffness).
`line.maxwell_branches`: `[{"K_N": ..., "tau_s": ...}, ...]`, each a spring K in series with a dashpot (relaxation time tau) in parallel with the curve: tension = curve(eps) + sum K (eps - alpha), d alpha/dt = (eps - alpha)/tau.
The relaxed (static) stiffness is the curve; the high-rate stiffness adds sum K. The static solve uses the relaxed state. Complex modulus of a branch: storage K w^2 tau^2 / (1 + w^2 tau^2), loss K w tau / (1 + w^2 tau^2).

### Instantaneous free-surface submergence
`waves.instantaneous_surface` (platform runner, default false): the submergence of line nodes and point elements (weight/buoyancy blend, drag, added mass) follows the instantaneous surface z = surface_z + eta(x, y, t) of the wave field instead of the still-water level. The static solve always uses the still-water level. In code: `Environment::elevation` (`WaveField::asElevation()`).

### Implicit internal damping
`numerics.implicit_damping` (single-line and platform runners, default false): treats the internal axial damping c_int implicitly (Strang splitting, see assumptions item 45). Removes the c_int time-step limit dt <= cfl m_l l0^2 / (2 c_int); only has an effect when `line.internal_damping_Ns` > 0.

### Non-uniform lines (sections)
`line.sections`: array of sections from the anchor to the fairlead, replacing `line.length_m` / `line.segments` (the line then has the sum of the sections' lengths and segments). Each section needs `length_m` and `segments` (segment length = length / segments) and may set `EA_N`, `mass_per_length_kg_m`, `density_kg_m3` or `weight_per_length_N_m`, `hydro_diameter_m`, `soil_diameter_m`, `nominal_area_m2`, `Cm`, `Cdt`, `Cdn`, `internal_damping_Ns`, `tension_curve`, `maxwell_branches`; keys that are omitted fall back to the line-level value (mass per length must be given somewhere). Gravity, soil (`soil.*`), the environment and the numerics stay global. Point elements can be placed by `arclength_m`, which now means the unstretched arc length from the anchor to the nearest node. Example: `examples/sections/chain_rope_chain.json`.
