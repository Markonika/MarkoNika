# Validation against Lopez-Olocco et al. (2022)

Source: Lopez-Olocco, T. et al., J. Mar. Sci. Eng. 10(5):676. A 27 m studless chain (EA 3.416e5 N, 0.0678 kg/m, 6.5 m depth, anchor 25 m from the
fairlead) with the fairlead driven on a horizontal circle (A = 0.125-0.225 m, T = 2.8-5.5 s); configurations WO_CW (no clump), CW1 (clump
122 g / 80 cm³ at Ls/3 from the fairlead), CW2 (at Ls/2). Reference data: Tables 8 (max) and 9 (min) fairlead tension, transcribed from the
paper's text layer into `examples/lopezolocco/measured.csv`. Configs: `examples/lopezolocco/{wocw,cw1,cw2}.json`; script `scripts/lopez_compare.py`.

## Results (60 segments clump-free, 81 with clump; max = mean of per-step maxima of the last 4 cycles)

| Config | Cases | Max tension, mean diff | mean abs diff | worst | Min tension, mean abs diff | worst |
|---|---|---|---|---|---|---|
| WO_CW | 35 (full grid) | +2.04 % | 2.04 % | +4.83 % (A 0.225, T 2.8) | 0.071 N | 0.215 N |
| CW1 | 9 (A 0.125/0.175/0.225 x T 2.8/3.5/5.5) | +0.76 % | 1.21 % | +2.32 % | 0.453 N | 0.911 N |
| CW2 | 9 (same) | -1.81 % | 1.81 % | -3.70 % | 0.174 N | 0.372 N |

Load-cell accuracy is quoted as about ±0.18 N, so the clump-free minima are within instrument accuracy; CW1 minima are over-predicted by
up to 0.9 N (the worst at A = 0.175, T = 2.8, short period, where the clump is near the touchdown region).

Static pretension at 25 m: WO 11.48 N, CW1 12.17 N, CW2 12.27 N. The paper states CW2 is about +10 % above WO; the model gives +6.9 %,
so the clump's static effect is somewhat **under**-predicted. The reason has not been isolated (candidates: clump volume/mass rounding,
chain axial stiffness, the soil/touchdown location).

## Assumptions and honest caveats
- Soil and damping parameters are the paper's Table 7 (Ksc 20, Dsc 0.1, friction 0); internal damping beta = 0.0007 s (c_int = beta·EA) as given there. No coefficient was tuned.
- The clump tangential drag coefficient (1.17) is **assumed** (not found in the paper). Sensitivity at CW1, A = 0.225, T = 2.8 (baseline 24.87 N,
  measured 24.31): Cd_t 0.8 -> 24.78 N; 1.5 -> 24.96 N; no added mass -> 25.00 N. The result is insensitive to these choices (≤0.5 %).
- Clump added mass rho·D³/3 is applied in the direction normal to the chain as in the paper's Eq. 25 (new anisotropic point element).
- The clump arclength is Ls/3 and Ls/2 with Ls = 14.33 m from the static solution, rounded to a node (81 segments); not converged in segment count for the clump cases.
- One paper value looks like a typo (A = 0.150, CW2, T = 5.0, max 15.72 N) and was not used in the subset runs.
- Only 53 of the 105 cases were run. Each took about 1-2 min with the explicit damping because the stability limit forces dt ~ 1e-5 s; with `numerics.implicit_damping` (assumptions item 45) the same results are obtained about 4-5 times faster (the validation test uses it).

Run: `ctest -R lopezolocco_validation` (four dynamic cases, skipped from the default suite); the default suite has a static test.
