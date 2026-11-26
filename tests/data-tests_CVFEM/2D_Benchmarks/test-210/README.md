<!--
SPDX-FileCopyrightText: © 2026 The openCSMP project

SPDX-License-Identifier: LGPL-3.0-only
-->

# Axisymmetric heat transfer (CVFEM_tests_main, 2D tests test-210a-c)

Axisymmetric models with `AxisymmetricThickness.h`: each element's out-of-plane thickness is factor · r.

## Conduction through an annulus (test-210a)

The cross-section 10 m ≤ r ≤ 210 m, 400 m high, of an annulus about the axis x = 0 (on the left boundary, factor 2π).
Inner wall 150 °C, outer wall 50 °C, top and base insulated, conduction only. From a uniform 50 °C start the
temperature approaches the steady radial-conduction profile T(r) = T₁ + (T₂ − T₁) ln(r/r₁)/ln(r₂/r₁); a plane model
would give a straight line instead, up to 34 K away. The test checks the revolved bulk volume against π(r₂² − r₁²)H
and the nodal temperature after 4000 yr against the logarithmic profile (limit 2 % of the wall difference; result
0.6 K). See `tests/CVFEM/CVFEM_AxisymmetricAnnulus_Test.h`.

## With radial outflow (test-210b)

The annulus 50 m ≤ r ≤ 250 m, inner wall 150 °C and 10.5 MPa, outer wall 50 °C and 10 MPa, so hot water flows
outward at about 3.15 kg/s. With the Péclet number A = c ṁ / (2π H λ) (here 2.64) the steady profile is
T(r) = T₁ + (T₂ − T₁)(r^A − r₁^A)/(r₂^A − r₁^A). The test measures ṁ from the steady pressure field and c from the
wall enthalpies. The plane conduction profile is 25 K away. Result after 1500 yr: 0.59 K (limit 2 K).

## A solid cylinder with the axis inside the mesh (test-210c)

The full section −200 m ≤ x ≤ 200 m of a solid cylinder about the axis x = 0, 400 m high (factor π: each half sweeps
180°). 30 of the 776 triangles straddle the axis; each part of such an element revolves about its own centroid. The
cylinder starts at 150 °C, its surface (both side walls) is held at 50 °C, top and base are insulated. The exact
temperature is the Bessel series T(r, t) = T₂ + (T₀ − T₂) Σ 2/(λₙ J₁(λₙ)) J₀(λₙ r/R) exp(−λₙ² Fo), λₙ the zeros of J₀.
After 250 yr (Fo = 0.14): revolved volume exact to 1.5 × 10⁻¹⁶ (with the barycentre weight alone 1.0 × 10⁻³),
temperature within 0.88 K (limit 2 K); the plane (slab) solution is 19.5 K away. See
`tests/CVFEM/CVFEM_AxisymmetricCylinder_Test.h`.

## References

- Conduction: Carslaw & Jaeger (1959), *Conduction of Heat in Solids*, 2nd ed., ch. 7 (§7.4 the solid cylinder,
  §3.3 the slab); Incropera et al. (2007), *Fundamentals of Heat and Mass Transfer*, 6th ed., §3.3.
- Radial flow: the profile follows by integrating the steady radial advection-conduction equation twice. It is the
  steady, finite-annulus counterpart of the radial thermal-injection problem of Avdonin (1964), a verification case
  for porous-flow codes (Ross et al. 1982, NUREG/CR-3097), and the radial analogue of the plane profile of
  Bredehoeft & Papadopulos (1965, *Water Resources Research* 1, 325–328).

## Meshes

- `test-210-unstructured`, `test-210-flow-unstructured`: the unstructured test-205 mesh (gmsh Frontal-Delaunay,
  target size 25 m) scaled by 0.4 and shifted by 10 m (test-210a) or 50 m (test-210b) in x.
- `test-210-cylinder`: gmsh Frontal-Delaunay triangles, target size 23 m, written by `make_test210_cylinder_mesh.py`.
