// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TWO_POINT_TRANSMISSIBILITY_H
#define TWO_POINT_TRANSMISSIBILITY_H

/** @file TwoPointTransmissibility.h
 *  @author Samuel Scott
 *  @brief Two-point (edge) transmissibilities of the CVFEM facets of a linear simplex.
 *  @details For the facet of edge (a, b) of a linear simplex (bar, triangle, tetrahedron) the element stiffness
 *  coefficient T_ab = -V_e grad N_a . K grad N_b (diagonal permeability tensor, as KgradP) gives the flux
 *  T_ab (Phi_a - Phi_b) from a to b, with Phi = p - rho g.x the phase potential. For a linear pressure the node totals
 *  equal those of the element-gradient (multi-point) facet fluxes. Unlike those, the flux between a and b depends on
 *  these two nodes only, so a phase in its own hydrostatic state exchanges nothing, whatever the other nodes of the
 *  element hold (an air flux next to a water table does not pick up the liquid's hydrostatic gradient). T_ab is
 *  positive on non-obtuse triangles and Delaunay tetrahedra with isotropic permeability, and can be negative
 *  otherwise. Per unit thickness in 2-D.
 *  @date 20.11.2026
 */

#include "Element.h"

#include <array>
#include <cmath>
#include <vector>

namespace csmp {

/// Two-point transmissibility of every facet of e, indexed as e.FV()->FacetEdgeNodes. Returns false and leaves
/// transmissibility unchanged for elements that are not linear simplices with one facet per edge.
template <uint32_t dim>
bool twoPointTransmissibilities(const Element<dim>& e, const TensorVariable<dim>& permeability,
                                std::vector<double>& transmissibility)
{
  constexpr uint32_t simplex_nodes{dim + 1U};
  constexpr uint32_t simplex_edges{dim * (dim + 1U) / 2U};
  if (e.Nodes() != simplex_nodes || e.FV()->Facets() != simplex_edges)
    return false;

  // grad N_k (k >= 1) = row k-1 of J^-1, grad N_0 = -sum; J[r][k] = x_{k+1}[r] - x_0[r]
  std::array<std::array<double, dim>, dim> m{};
  for (uint32_t k{0U}; k < dim; k++)
    for (uint32_t r{0U}; r < dim; r++)
      m[r][k] = e.N(k + 1U)->Coordinate()[r] - e.N(0U)->Coordinate()[r];
  std::array<std::array<double, dim>, dim> inverse{};
  double determinant{0.};
  if constexpr (dim == 1U) {
    determinant = m[0][0];
    inverse[0][0] = 1. / determinant;
  }
  else if constexpr (dim == 2U) {
    determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    inverse[0][0] = m[1][1] / determinant;
    inverse[0][1] = -m[0][1] / determinant;
    inverse[1][0] = -m[1][0] / determinant;
    inverse[1][1] = m[0][0] / determinant;
  }
  else {
    determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                  + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    inverse[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / determinant;
    inverse[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / determinant;
    inverse[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / determinant;
    inverse[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / determinant;
    inverse[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / determinant;
    inverse[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / determinant;
    inverse[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / determinant;
    inverse[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / determinant;
    inverse[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / determinant;
  }
  std::array<std::array<double, dim>, simplex_nodes> gradient{};
  for (uint32_t k{1U}; k < simplex_nodes; k++)
    for (uint32_t r{0U}; r < dim; r++) {
      gradient[k][r] = inverse[k - 1U][r];
      gradient[0][r] -= gradient[k][r];
    }
  const double volume{std::fabs(determinant) / (dim == 3U ? 6. : (dim == 2U ? 2. : 1.))};

  transmissibility.assign(simplex_edges, 0.);
  for (uint32_t i{0U}; i < simplex_edges; i++) {
    uint32_t a{0U}, b{0U};
    e.FV()->FacetEdgeNodes(i, a, b);
    for (uint32_t r{0U}; r < dim; r++)
      transmissibility[i] -= volume * gradient[a][r] * permeability(r, r) * gradient[b][r];
    // Negative on the edge opposite an obtuse angle: T (Phi_a - Phi_b) would then run from low to high potential, and
    // with the air weighting of the pressure equation it fills a saturated node below the water table with phantom air
    // volume until it drains empty. No two-point flux on such an edge; the other edges keep the nodes coupled.
    if (transmissibility[i] < 0.) transmissibility[i] = 0.;
  }
  return true;
}

} // namespace csmp

#endif // TWO_POINT_TRANSMISSIBILITY_H
