# Implications of the FOWT multiphysics-tools review for this code

Source: M. Nika, *Multiphysics Modelling Tools for Floating Offshore Wind Turbines: A Systematic Literature Review* (preprint, `FOWT_Review_Multiphysics_Tools.tex`).
Statements about other codes below are **as stated in that review**; nothing here was tested against, or copied from, MoorDyn, OpenFAST or any other code.

## Where this code sits
In the review's taxonomy this is a **Tier VI lumped-mass dynamic mooring code** (the class of MoorDyn v2 and MoodyMarine), with a coupled 6-DOF body and
frequency-dependent hydrodynamics (Cummins) in addition, i.e. partly Tier I-like at the platform level. The review's fidelity guidance (Sec. "Fidelity trade-off"): lumped-mass is
the recommended level for wave-frequency response, fatigue, synthetic-rope material nonlinearity and shared moorings; nonlinear FEM only for deep water (> 500 m), torsion/bending-critical
cases and certification. That matches the scope of the code.

## Feature gap list (review-driven), against the current code
| Review topic | Review statement | This code today | Priority |
|---|---|---|---|
| Large platform rotations | wave-frequency surge-pitch coupling dominates semi-submersibles; OC5: pitch under-predicted in steep waves | linear body dynamics (constant A, B, C, omega = d theta/dt) with exact rotation of the fairlead positions; nonlinear hydrostatics and 3-D rotational dynamics are not modelled. The CON2 excess is **not** caused by this (see `bem_paredes.md`, quasi-static pitch curve) | 3 (for 3-D motions and large angles of other platforms; no longer tied to CON2) |
| Second-order wave loads | difference-frequency (slow drift) mandatory for spar, recommended for semi; sum-frequency for TLP; OC5/OC6 main open gap; OpenFAST gained Newman approximation and full QTF import | first-order only (`wave_force_table`, Cummins); no mean drift, no QTF | **2** (Newman approximation first, QTF table import second) |
| Viscous body damping | OC6 Ph. I: viscous drag on heave plates adds 15-25 % heave damping; potential-flow-only models miss it | body quadratic damping `Dq` exists but is zero in all validation cases (the Paredes thesis gives no values); radiation B only | 3 (needs data; do not guess coefficients) |
| Shared mooring / multi-body | lines between adjacent floaters (about 30 % fewer anchors) are supported by MoorDyn v2 and are an emerging design trend | one rigid body per case | 3 |
| One-line-broken ALS cases | IEC 61400-3-2 requires them for semi-submersibles | no line removal during a run | 4 (cheap: remove a line's force at time t) |
| 6-DOF buoys / clumps | MoorDyn v2 and OrcaFlex treat clumps and buoys as rigid bodies | point masses with drag (our CON2/CON1 point elements) | 4 |
| Seabed | MoorDyn v2: bathymetry grid and Coulomb friction | flat bed, spring-damper plus Coulomb friction | 4 (bathymetry only) |
| Fibre-rope material | viscoelastic model "essential" above 300-500 m | implemented (static curve + Maxwell branches), **never validated against rope data** | validation, see below |
| Bending / torsion | RIFLEX/OrcaFlex: bending, torsion; MoorDyn v2: bending stiffness | bending yes (no torsion) | 5 (torsion only for dynamic cables) |
| Nonlinear Froude-Krylov | MoodyMarine: improved accuracy in steep waves | Airy + Wheeler, Morison | 5 |
| Fatigue post-processing | OrcaFlex offers S-N and fracture-mechanics fatigue; mooring fatigue is a central driver | none | 3 (rainflow counting + S-N damage on tension series is cheap) |
| Quasi-static screening | two-stage workflow (QS sweep, then dynamic DLCs) saves 60-80 % of simulation time | static solves exist (`relaxStatic`, catenary start), no screening interface | 4 (restoring curve / stiffness sweep tool) |
| Coupling interfaces | MoorDyn-C is used as a library for OpenFOAM and others; FMI 3.0 is the standard co-simulation interface; loose coupling needs dt <~ 0.01 s | explicit and implicit coupling inside the code; no external API | 3 (C API + FMU wrapper, if the code is to be used with other tools) |
| Wind / aerodynamic loads | AHSE tools couple rotor loads | external body force only (`F0`, tables) | outside scope unless coupled |

## Validation targets named by the review
* OC4 DeepCwind semi-submersible (code-to-code reference geometry), OC5 Phase II (MARIN 1:50 model tests, the main FOWT code-to-experiment dataset), OC6 Phase II VolturnUS-S (UMaine 1:70, with experimental uncertainty bounds), OC6 Phase IV TetraSpar (1:43),
  the SINTEF WINDMOOR chain-fibre-chain case (QS peak tension off by more than a factor 2 against RIFLEX FEM).
* These are the datasets that would give the rope law and multi-line platforms their first test. They need the data hosts listed in earlier notes (`mhkdr.openei.org`, `data.openei.org`, 4TU/TU Delft hosts) to be allowed; published geometry/mooring parameters may be used, source files of other codes must not be copied.
* OC6 Phase II introduced experimental **uncertainty bounds**; our Paredes/Azcona comparisons use by-eye values and should report uncertainty the same way (already flagged in the paper).

## Suggested order
1. Finish the CON2 investigation: linear dynamic stiffness of the legs at the wave frequency, then drag/inertia switched off one at a time (the large-angle body was tested as an explanation and refuted).
2. Newman-approximation slow-drift force from a mean-drift table (Capytaine can supply first-order data; mean-drift/QTF data must come from a BEM that provides it), then full QTF table import.
3. Fatigue post-processing (rainflow, S-N), line-removal (ALS) option, C API / FMU wrapper.
4. Multi-body and shared mooring; 6-DOF buoys; bathymetry; quasi-static screening tool.
5. Validation against OC5/OC6 data once the hosts are reachable (rope law, multi-line platforms).

## Positioning in the paper
The paper can cite the review for the taxonomy and add a feature-comparison table against MoorDyn v2/MoodyMarine/OrcaFlex **based on the review's descriptions only**
(not on tests); it must keep stating that no code-to-code comparison was made.
