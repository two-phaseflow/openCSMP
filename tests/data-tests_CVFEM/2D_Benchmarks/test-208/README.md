<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Isothermal hydrostatic air column (CVFEM_tests_main, 2D tests test-208a-b)

An air column 1000 m × 2000 m at 20 °C (top 1 atm fixed, no-flow base and sides, permeability 10⁻¹³ m²), with the
pore water at the 2 % floor that the H2O-NaCl equation of state always keeps and the rest of the pore space filled
with air. The exact solution is no flow, with the exponential pressure profile p(z) = p_top exp(−g (z − z_top) / (R_a T))
of the code's ideal-gas air; the air density falls by 21 % from base to top. See
`tests/CVFEM/CVFEM_GasHydrostaticColumn_Test.h` for the set-up and the pass criteria.

The meshes are the test-205 meshes scaled by 2 (test-208a structured, test-208b unstructured; facet-mean density in
the gravity term):

- `test-208-structured`: 20 × 40 squares split into triangles with alternating diagonals (50 m).
- `test-208-unstructured`: gmsh Frontal-Delaunay triangles, target size 50 m.

test-209a-b (closed-box air compression, `tests/CVFEM/CVFEM_AirCompressionBox_Test.h`) use the test-207 mesh.
