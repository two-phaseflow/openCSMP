// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_RECHARGE_MOUND_TEST_H
#define CVFEM_RECHARGE_MOUND_TEST_H

/** @file CVFEM_RechargeMound_Test.h
 *  @author Samuel Scott
 *  @brief Rain, vadose zone and water table of the CVFEM PHX scheme: the steady groundwater mound under uniform
 *  recharge against the Dupuit solution.
 *  @details Isothermal box 0 <= x <= L, base at y = -D, land surface at y = 0, air in the hydraulic-mixture mode with
 *  rain (CVFEM_PHX_Scheme::ActivateAir(0., false, true, rain), OpenBoundaries). The land surface (TOP) is held at
 *  atmospheric pressure and receives the rain. The RIGHT side is a fixed water level h0 above the base: the test holds
 *  its nodes at p = p_atm + rho_l g max(0, h0 - (y + D)), hydrostatic below the water level and atmospheric above it
 *  (a seepage face). LEFT (the groundwater divide) and the base are closed. The box starts with the water table flat
 *  at the outlet water level, water-saturated and hydrostatic below it, and at the residual liquid saturation (0.3)
 *  with air at atmospheric pressure above it (CVFEM_AirPhaseLoader); the rain then builds the mound. Starting
 *  water-saturated and draining instead makes the air time-step limit collapse (air entering nodes that hold almost
 *  none moves at a very high pore velocity), which is too slow for the test suite. With the Dupuit approximation the steady water
 *  table height above the base is
 *      h(x)^2 = h0^2 + (W / K) (L^2 - x^2),
 *  with W the recharge (volume of groundwater per unit area and time) and K = k rho_l g / mu_l the hydraulic
 *  conductivity. The test takes W from the rain mass the scheme actually assigns (sum of "rain recharge top" times
 *  "pore volume", divided by L and rho_l) and K from the liquid density and viscosity at the base. The model water
 *  table in each mesh column is where the fluid pressure crosses p_atm, interpolated linearly between nodes.
 *
 *  With rain 0 (test-211b) the water table stays at the outlet level: the exact solution is no flow, and the same
 *  criteria check that the water table has not moved.
 *
 *  Pass criteria, after the configured duration:
 *  - the water table at the divide (x = 0) within max_divide_error_ of the Dupuit height;
 *  - the water table in every column within max_column_error_ of the Dupuit height, except the column at the
 *    outlet, where the seepage face lifts the water table above the Dupuit curve;
 *  - for information: the change of the water table at the divide over the last 10 % of the run (steady state).
 *  It writes x, the model and the Dupuit water table per column to a CSV file.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_2D_VVCase and tests/CVFEM/CVFEM_HydrostaticColumn_Test.
 *
 *  Reference: the Dupuit-Forchheimer approximation for unconfined flow with uniform recharge, e.g. Fetter C.W.
 *  (2001), Applied Hydrogeology, 4th ed., Prentice Hall, section 4.14; Bear J. (1972), Dynamics of Fluids in Porous
 *  Media, Elsevier.
 *  @date 20.11.2026
 */

#include "Test.h"

#include "ComputationalSettings.h"
#include "CVFEM_PHX_Scheme.h"
#include "PermeabilityVisitor.h"
#include "CVFEM_AirPhaseLoader.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace csmp {

/// Steady groundwater mound under uniform rain against the Dupuit solution.
class CVFEM_RechargeMound_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt); rain in mm/yr,
    /// outlet_water_level h0 in m above the base.
    CVFEM_RechargeMound_Test(const std::string& mesh_name, const std::string& config_name, double rain_mm_per_year,
                             double outlet_water_level);
    ~CVFEM_RechargeMound_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, builds the CVFEM PHX scheme with air and rain, permeability tensor and pore volumes.
    void CreateModel();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Pressure hydrostatic below the outlet water level and atmospheric above it on all nodes, with the RIGHT nodes
    /// held at it; liquid density from the equation of state (two passes).
    void InitializeFlatWaterTable();
    /// Air above the water table: liquid saturation vadose_saturation_, air at the local pressure in the rest.
    void LoadAirAboveWaterTable();
    /// Steps to max_time.
    void Advance(double max_time);
    /// Model water table height above the base per mesh column (x -> h), from the p_atm crossing.
    std::map<double, double> WaterTable() const;
    /// Rain mass rate assigned by the scheme [kg/s per m out of plane].
    double RainMassRate() const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    double rain_mm_per_year_;
    double outlet_water_level_;
    uint64_t steps_;
    double x_min_;
    double x_max_;
    double y_base_;
    double y_top_;
    double max_divide_error_;
    double max_column_error_;
    double vadose_saturation_;
    double air_weight_;   ///< rho_air(T, p_atm) g [Pa/m], air at rest

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    std::unique_ptr<CVFEM_AirPhaseLoader> air_phase_loader_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr double gravity_{9.80665};
    static constexpr double atmospheric_pressure_{101325.};
};

} // namespace csmp

#endif // CVFEM_RECHARGE_MOUND_TEST_H
