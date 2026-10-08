# Assumptions and open items

1. **K = EA.** Paredes (2016) Eqs. 3.54-3.55 contain `m_l g s / K`, which equals `H s / EA`, so the
   "stiffness K = 200 kN/m" of the sec. 3.5.1 static test is the axial stiffness EA = 200 kN.
   The test additionally runs EA = K*L (stiff) as a diagnostic.
2. Static test is in air: weight per length w = m_l * g (dry).
3. Static relaxation uses fictitious nodal mass (massFactor * EA / l0, pseudo-dt = 1); only the
   converged equilibrium is physical. massFactor = 2 was unstable, 32 is the default.
4. Reported L2 error is relative, over nodes at Lagrangian coordinate s_i = i*l0.
5. Still needed: fairlead coordinates and line data for the Paredes buoy (thesis Figs 5.4-5.5,
   Tables 5.3, 5.5-5.10) - to be read from the supplied thesis PDF at milestone 7.
6. `chalmers_model_inputs.json` and `chalmers_table7_max_tension.csv` have not been supplied yet.
