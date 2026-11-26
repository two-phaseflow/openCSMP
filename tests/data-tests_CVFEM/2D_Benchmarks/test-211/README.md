<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Water table with air above it (CVFEM_tests_main, 2D tests test-211a-c)

A box 600 m × 200 m, isothermal (20 °C), air in the hydraulic-mixture mode with rain on the land surface (TOP, 1 atm).
The RIGHT side is a fixed water level 100 m above the base: hydrostatic water below it, air at rest above it (a seepage
face). The divide (LEFT) and the base are closed. The test starts with the water table flat at the outlet level, water
at rest below it and air at rest above it (residual liquid saturation 0.3). See `tests/CVFEM/CVFEM_RechargeMound_Test.h`.

- **test-211a, static water table:** no rain, 0.5 yr. Exact solution: no flow. Result: water table within
  8 × 10⁻⁵ (relative) of 100 m at the divide, 1.4 × 10⁻⁴ in every column.
- **test-211b, recharge mound:** rain 0.96 m/yr, permeability 2 × 10⁻¹³ m², 25 yr. The steady water table is compared
  with the Dupuit solution h² = h₀² + (W/K)(L² − x²), about 125 m at the divide. Result: 0.75 % at the divide, 1.9 % in
  the worst column (limits 2 % and 5 %).

- **test-211c, static water table on a mesh with obtuse triangles:** as test-211a on `test-211-obtuse`, the same box
  with the node row 10 m below the water table shifted sideways by ±8 m, alternately. 29 triangle edges with a
  negative finite-element coefficient join a water-table node to a node below it. Same criteria as test-211a.

All three need the two-point air fluxes (`TwoPointTransmissibility.h`, `CVFEM_TwoPointUpwind_Operators.h`): with the element
pressure gradient the air in a triangle that straddles the water table picks up the liquid's hydrostatic gradient and the
water table is not at rest. On the edge opposite an obtuse angle the two-point coefficient is negative; such edges
carry no two-point flux (`TwoPointTransmissibility.h`), otherwise the air-weighted pressure equation pumps phantom air
into the saturated node below the water table and drains it (test-211c).

Reference: the Dupuit-Forchheimer approximation with uniform recharge, Fetter (2001), *Applied Hydrogeology*, 4th ed.,
§4.14; Bear (1972), *Dynamics of Fluids in Porous Media*.

The mesh `test-211-structured` is 30 × 20 rectangles (20 m × 10 m) split into triangles with alternating diagonals;
`make_test211_obtuse_mesh.py` writes `test-211-obtuse` (with a zero shift it reproduces `test-211-structured`).
