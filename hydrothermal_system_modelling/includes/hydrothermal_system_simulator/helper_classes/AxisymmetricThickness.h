// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef AXISYMMETRIC_THICKNESS_H
#define AXISYMMETRIC_THICKNESS_H

/** @file AxisymmetricThickness.h
 *  @author Samuel Scott
 *  @brief Axisymmetric 2-D models: element "thickness" = 2 pi r.
 *  @details The CVFEM operators (capacitances, conductances, gravity terms, finite-volume transport) and the
 *  PoreVolumeVisitor all multiply by the element property "thickness", the out-of-plane extent of a 2-D model.
 *  Revolving the (x, y) half-plane about the vertical axis x = axis_x turns that extent into 2 pi r, so
 *      thickness_e = factor |x_barycentre - axis_x|,   factor = 2 pi for a full revolution,
 *  gives every element its exact revolved volume (Pappus: the integral of 2 pi r dA is 2 pi r_centroid A). The factor
 *  is uniform per element, so the CVFEM control volumes get equal thirds of each revolved element volume (a
 *  second-order approximation of the true revolved sub-volumes).
 *
 *  Linear triangles only. Elements cut by the axis are split along it (as CrustalFate's HalfRadialVolumeVisitor);
 *  with the axis inside the mesh use factor pi (each half sweeps 180 degrees). Boundary sources that are not integrated with "thickness" (nodal heat or mass fluxes set per boundary
 *  face) must be weighted with RadialWeight() at the face.
 *  @date 20.11.2026
 */

#include "Model.h"
#include "Element.h"

#include <array>
#include <cmath>
#include <vector>

namespace csmp {

/// Axis and revolution factor of an axisymmetric 2-D model.
struct AxisymmetricSettings {
  bool enabled{false};
  double axis_x{0.};          ///< x of the symmetry axis [m]
  /// Revolution angle: 2 pi for a half-section with the axis on its boundary; pi for a full section with the axis
  /// inside it (each half sweeps 180 degrees, together one revolution: CrustalFate's HalfRadialVolumeVisitor,
  /// radial_factor 1)
  double factor{2. * M_PI};

  /// Out-of-plane extent at x: factor |x - axis_x| if enabled, else 1.
  double RadialWeight(double x) const;
};

inline double AxisymmetricSettings::RadialWeight(double x) const
{
  return enabled ? factor * std::fabs(x - axis_x) : 1.;
}

/// Integral of |x - axis_x| over a triangle; for a triangle cut by the axis the two parts are integrated separately
/// (each part's area times the distance of its own centroid), else -1 (not cut: the barycentre weight is exact).
inline double revolvedMoment(const std::vector<std::array<double, 2>>& corners, double axis_x)
{
  bool left{false}, right{false};
  for (const auto& c : corners) { left |= c[0] < axis_x - 1.e-9; right |= c[0] > axis_x + 1.e-9; }
  if (!(left && right)) return -1.;
  double moment{0.};
  for (const double side : {-1., 1.}) {
    std::vector<std::array<double, 2>> part;   // the triangle clipped to one side of the axis (Sutherland-Hodgman)
    for (size_t i{0U}; i < corners.size(); ++i) {
      const auto& a{corners[i]}; const auto& b{corners[(i + 1) % corners.size()]};
      const bool a_in{side * (a[0] - axis_x) >= 0.}, b_in{side * (b[0] - axis_x) >= 0.};
      if (a_in) part.push_back(a);
      if (a_in != b_in) {
        const double t{(axis_x - a[0]) / (b[0] - a[0])};
        part.push_back({axis_x, a[1] + t * (b[1] - a[1])});
      }
    }
    double area{0.}, cx{0.};
    for (size_t i{0U}; i < part.size(); ++i) {
      const auto& a{part[i]}; const auto& b{part[(i + 1) % part.size()]};
      const double cross{a[0] * b[1] - b[0] * a[1]};
      area += cross; cx += (a[0] + b[0]) * cross;
    }
    if (std::fabs(area) > 0.) moment += std::fabs(0.5 * area) * std::fabs(cx / (3. * area) - axis_x);
  }
  return moment;
}

/// Stores "thickness" = settings.RadialWeight(x_barycentre) on every element of region "Model" (2-D only) and
/// returns the total revolved bulk volume, for a check against the analytic volume of the domain.
template <uint32_t dim>
double setAxisymmetricThickness(Model<dim>& model, const AxisymmetricSettings& settings)
{
  static_assert(dim == 2U, "setAxisymmetricThickness: 2-D models only");
  const Index thickness_key{model.Database().StorageKey("thickness")};
  double total_volume{0.};
  const Region<dim>& region{model.Region("Model")};
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it) {
    Element<dim>& element{**element_it};
    double weight{settings.RadialWeight(element.BaryCenter()[0])};
    if (settings.enabled && element.Nodes() == 3U) {
      // an element cut by the axis: each part revolves about its own centroid (CrustalFate HalfRadialVolumeVisitor)
      std::vector<std::array<double, 2>> corners;
      for (uint32_t i{0U}; i < 3U; ++i) corners.push_back({element.N(i)->x(), element.N(i)->y()});
      const double moment{revolvedMoment(corners, settings.axis_x)};
      if (moment >= 0. && element.Volume() > 0.) weight = settings.factor * moment / element.Volume();
    }
    element.Store(thickness_key, ScalarVariable(ANY, weight));
    total_volume += weight * element.Volume();
  }
  return total_volume;
}

} // namespace csmp

#endif // AXISYMMETRIC_THICKNESS_H
