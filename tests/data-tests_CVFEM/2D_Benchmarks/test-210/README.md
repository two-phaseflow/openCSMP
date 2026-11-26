<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Steady heat transfer in an axisymmetric annulus (CVFEM_tests_main, 2D tests test-210a-d)

The cross-section 10 m ≤ r ≤ 210 m, 400 m high, of an annulus about the axis x = 0, made axisymmetric with
`AxisymmetricThickness.h` (element thickness 2π r). Inner wall 150 °C, outer wall 50 °C, top and base insulated,
conduction only. From a uniform 50 °C start the temperature approaches the steady radial-conduction profile
T(r) = T₁ + (T₂ − T₁) ln(r/r₁)/ln(r₂/r₁); a plane model would give a straight line instead, up to 34 K away.
The test checks the revolved bulk volume against π(r₂² − r₁²)H and the nodal temperature after 4000 yr against
the logarithmic profile (limit 2 % of the wall difference; the result is 1.0 K on the structured and 0.6 K on the unstructured mesh, mostly next to the inner wall). See `tests/CVFEM/CVFEM_AxisymmetricConduction_Test.h`.

## With radial outflow (test-210c/d)

The annulus 50 m ≤ r ≤ 250 m, inner wall 150 °C and 10.5 MPa, outer wall 50 °C and 10 MPa, so hot water flows
outward at about 3.15 kg/s. With the Péclet number A = c ṁ / (2π H λ) (here 2.64) the steady profile is
T(r) = T₁ + (T₂ − T₁)(r^A − r₁^A)/(r₂^A − r₁^A). The test measures ṁ from the steady pressure field and c from the
wall enthalpies. The plane conduction profile is 25 K away, the axisymmetric conduction profile further still.
Result after 1500 yr: 0.73 K on the structured and 0.59 K on the unstructured mesh (limit 2 K).

## References

- Conduction: Carslaw & Jaeger (1959), *Conduction of Heat in Solids*, 2nd ed., ch. 7; Incropera et al. (2007),
  *Fundamentals of Heat and Mass Transfer*, 6th ed., §3.3.
- Radial flow: the profile follows by integrating the steady radial advection-conduction equation twice. It is the
  steady, finite-annulus counterpart of the radial thermal-injection problem of Avdonin (1964), a verification case
  for porous-flow codes (Ross et al. 1982, NUREG/CR-3097), and the radial analogue of the plane profile of
  Bredehoeft & Papadopulos (1965, *Water Resources Research* 1, 325–328).

## Meshes

The test-205 meshes scaled by 0.4 and shifted by 10 m (test-210a/b) or 50 m (test-210c/d) in x:

- `test-210-structured`, `test-210-flow-structured`: 20 × 40 squares split into triangles with alternating
  diagonals (10 m).
- `test-210-unstructured`, `test-210-flow-unstructured`: gmsh Frontal-Delaunay triangles, target size 10 m.
