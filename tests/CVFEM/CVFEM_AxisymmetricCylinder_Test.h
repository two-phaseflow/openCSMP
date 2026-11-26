// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_AXISYMMETRIC_CYLINDER_TEST_H
#define CVFEM_AXISYMMETRIC_CYLINDER_TEST_H

/** @file CVFEM_AxisymmetricCylinder_Test.h
 *  @author Samuel Scott
 *  @brief Axisymmetric 2-D geometry with the axis inside the mesh (AxisymmetricThickness.h, factor pi): a solid
 *  cylinder cooling by conduction.
 *  @details The mesh is the full section -R <= x <= R of a solid cylinder about the axis x = 0, with elements that
 *  straddle the axis. Both side walls (the cylinder surface r = R) are held at T2, top and base are closed and
 *  insulated, and the cylinder starts at T0. Conduction only. The exact temperature is the Bessel series
 *      T(r, t) = T2 + (T0 - T2) sum_n 2 / (l_n J1(l_n)) J0(l_n r / R) exp(-l_n^2 Fo),   Fo = kappa t / R^2,
 *  with l_n the zeros of J0 and kappa = lambda / C the bulk thermal diffusivity (C from the rock and the pore water
 *  at the end of the run). The test checks (i) the revolved bulk volume against pi R^2 H, which needs the elements
 *  cut by the axis split along it, and (ii) the nodal temperature after the run against T(r, t), within 2 % of
 *  T0 - T2. It reports the largest difference of the plane (slab) solution for comparison and writes x, y, T, the
 *  exact and the plane temperature per node to a CSV file.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_AxisymmetricAnnulus_Test.
 *
 *  Reference: Carslaw H.S. & Jaeger J.C. (1959), Conduction of Heat in Solids, 2nd ed., Oxford University Press,
 *  section 7.4 (the infinite cylinder with constant surface temperature) and section 3.3 (the slab).
 *  @date 26.11.2026
 */

#include "Test.h"

#include "AxisymmetricThickness.h"
#include "ComputationalSettings.h"
#include "CVFEM_PHX_Scheme.h"
#include "PermeabilityVisitor.h"

#include <memory>
#include <string>

namespace csmp {

/// A solid cylinder cooling by conduction, full section with the axis inside the mesh, against the Bessel series.
class CVFEM_AxisymmetricCylinder_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt).
    CVFEM_AxisymmetricCylinder_Test(const std::string& mesh_name, const std::string& config_name);
    ~CVFEM_AxisymmetricCylinder_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, sets the axisymmetric thickness (factor pi), builds the CVFEM PHX scheme,
    /// permeability tensor and pore volumes. Returns the total revolved bulk volume.
    double CreateModel();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Reads R, the height, the start temperature T0 (largest nodal value) and the wall temperature T2.
    void ReadGeometryAndTemperatures();
    /// Steps to max_time.
    void Advance(double max_time);
    /// Bulk thermal diffusivity lambda / C [m2/s], C averaged over the nodes.
    double Diffusivity() const;
    /// Exact temperature at radius r for the Fourier number fourier.
    double ExactTemperature(double r, double fourier) const;
    /// Plane (slab of half-width R) temperature at x for the Fourier number fourier.
    double PlaneTemperature(double x, double fourier) const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    uint64_t steps_;
    double radius_;
    double height_;
    double start_temperature_;
    double wall_temperature_;

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    ComputationalSettings run_settings_;
    AxisymmetricSettings axisymmetric_;

    static constexpr double year_{31536000.};
};

} // namespace csmp

#endif // CVFEM_AXISYMMETRIC_CYLINDER_TEST_H
