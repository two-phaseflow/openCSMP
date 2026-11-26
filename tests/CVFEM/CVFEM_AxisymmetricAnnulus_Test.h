// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_AXISYMMETRIC_ANNULUS_TEST_H
#define CVFEM_AXISYMMETRIC_ANNULUS_TEST_H

/** @file CVFEM_AxisymmetricAnnulus_Test.h
 *  @author Samuel Scott
 *  @brief Axisymmetric 2-D geometry (AxisymmetricThickness.h): steady heat transfer through a thick cylindrical
 *  shell, by conduction alone or with radial fluid flow.
 *  @details The mesh is the cross-section r1 <= x <= r2 of an annulus about the axis x = 0. The inner wall (LEFT) is
 *  held at T1 and p1, the outer wall (RIGHT) at T2 and p2, top and base are closed and insulated. Water flows
 *  radially at the mass rate m (kg/s through the whole shell), constant in r at steady state. With the Peclet number
 *      A = c m / (2 pi H lambda)
 *  (c the fluid heat capacity, H the height, lambda the bulk thermal conductivity), the steady energy balance
 *  c m dT/dr = 2 pi H lambda d/dr (r dT/dr) gives
 *      T(r) = T1 + (T2 - T1) (r^A - r1^A) / (r2^A - r1^A),
 *  which becomes the conduction profile T1 + (T2 - T1) ln(r / r1) / ln(r2 / r1) for p1 = p2 (A = 0). A plane model
 *  gives other profiles (a straight line without flow). The test measures m from the steady pressure field (Darcy
 *  flux of every element, r rho q_r averaged over the cross-section) and c from the liquid enthalpies of the two
 *  walls, (h1 - h2) / (T1 - T2). It checks (i) the total revolved bulk volume against pi (r2^2 - r1^2) H and
 *  (ii) the nodal temperature after the run against T(r), within 2 % of the wall difference (the element-centroid
 *  thickness integrates 1/r with the midpoint rule). It writes x, y, T, the exact profile, the conduction-only
 *  profile and the plane conduction profile per node to a CSV file.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_2D_VVCase.
 *
 *  References:
 *  - conduction (test-210a/b): steady radial conduction through a hollow cylinder, Carslaw H.S. & Jaeger J.C.
 *    (1959), Conduction of Heat in Solids, 2nd ed., Oxford University Press, chapter 7; Incropera F.P.,
 *    DeWitt D.P., Bergman T.L. & Lavine A.S. (2007), Fundamentals of Heat and Mass Transfer, 6th ed., Wiley,
 *    section 3.3 (radial systems);
 *  - radial flow (test-210c/d): T(r) above follows by integrating the steady radial advection-conduction equation
 *    twice. It is the steady, finite-annulus counterpart of the radial thermal-injection problem of Avdonin N.A.
 *    (1964), "Some formulas for calculating the temperature field of a stratum subject to thermal injection",
 *    Neft' i Gaz, a standard verification case for porous-flow codes (Ross B., Mercer J.W., Thomas S.D. & Lester B.H.
 *    (1982), Benchmark problems for repository siting models, NUREG/CR-3097), and the radial analogue of the plane
 *    profile of Bredehoeft J.D. & Papadopulos I.S. (1965), "Rates of vertical groundwater movement estimated from the
 *    Earth's thermal profile", Water Resources Research 1(2), 325-328.
 *
 *  @date 20.11.2026
 */

#include "Test.h"

#include "AxisymmetricThickness.h"
#include "ComputationalSettings.h"
#include "CVFEM_PHX_Scheme.h"
#include "PermeabilityVisitor.h"

#include <array>
#include <memory>
#include <string>

namespace csmp {

/// Steady radial heat transfer in an axisymmetric annulus, with or without radial flow, against the exact profile.
class CVFEM_AxisymmetricAnnulus_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt).
    CVFEM_AxisymmetricAnnulus_Test(const std::string& mesh_name, const std::string& config_name);
    ~CVFEM_AxisymmetricAnnulus_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, sets the axisymmetric thickness, builds the CVFEM PHX scheme, permeability
    /// tensor and pore volumes. Returns the total revolved bulk volume.
    double CreateModel();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Reads r1, r2, the height and the wall temperatures T1 (LEFT) and T2 (RIGHT) from the mesh and configuration.
    void ReadGeometryAndWallTemperatures();
    /// Steps to max_time.
    void Advance(double max_time);
    /// Radial mass rate through the shell [kg/s] (element Darcy fluxes, r rho q_r averaged over the cross-section)
    /// [0] and the largest relative deviation of an element's r rho q_r from that average [1].
    std::array<double, 2> RadialMassRate() const;
    /// Mean liquid enthalpy of the nodes on the wall at x.
    double WallEnthalpy(double x) const;
    /// Steady temperature at radius r for the Peclet number peclet_.
    double ExactTemperature(double r) const;
    /// Steady conduction-only temperature at radius r (ln r profile).
    double ConductionTemperature(double r) const;
    /// Plane (straight-line) conduction temperature at x.
    double PlaneTemperature(double x) const;
    /// Largest |T - ExactTemperature| over all nodes; writes the per-node CSV.
    double MaxTemperatureError() const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    uint64_t steps_;
    double inner_radius_;
    double outer_radius_;
    double height_;
    double inner_temperature_;
    double outer_temperature_;
    double peclet_;

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    ComputationalSettings run_settings_;
    AxisymmetricSettings axisymmetric_;

    static constexpr double year_{31536000.};
};

} // namespace csmp

#endif // CVFEM_AXISYMMETRIC_ANNULUS_TEST_H
