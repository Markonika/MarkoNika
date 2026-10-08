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
| `line.internal_damping_Ns` | c_int in T = EA eps + c_int d(eps)/dt (0 = off; also the slack regularisation, logged) |
| `environment.hydro`, `.seabed` | enable Morison forces + submerged weight; enable seabed |
| `environment.rho_w_kg_m3`, `.water_surface_z_m`, `.seabed_z_m` | water density, still-water level, flat seabed height |
| `soil.stiffness_Pa_per_m`, `.damping_factor`, `.friction`, `.v_lim_m_s` | Ks, zeta_s, mu_s, v_lim (Eqs. 3.32-3.36) |
| `anchor_m`, `fairlead_rest_m` | anchor (fixed) and top end at rest |
| `initial.shape` | `touchdown` (catenary lying on the bed, then static relaxation) |
| `motion.type` | `none` or `circle_xz` (top end on a circle about `centre_m`, default = rest position) |
| `motion.radius_m`, `.period_s`, `.direction` (+1 counter-clockwise in x-z, -1 clockwise), `.phase_deg`, `.ramp_cycles`, `.cycles` | excitation |
| `numerics.scheme` (`rk4`/`verlet`), `.cfl`, `.dt_s` (override), `.relax_force_tol`, `.relax_max_steps` | integration and static relaxation |
| `statistics.first_cycle` | first cycle included in the mean of the cycle maxima |
| `output.directory`, `.tag`, `.dt_out_s`, `.write` | output control |
