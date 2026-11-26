// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_GAS_HYDROSTATIC_COLUMN_TEST_H
#define CVFEM_GAS_HYDROSTATIC_COLUMN_TEST_H

/** @file CVFEM_GasHydrostaticColumn_Test.h
 *  @author Samuel Scott
 *  @brief Well-balancedness of the air phase of the CVFEM PHX scheme: an isothermal air column must stay at rest in
 *  its exact (exponential) hydrostatic profile.
 *  @details Air-phase counterpart of tests/CVFEM/CVFEM_HydrostaticColumn_Test (test-205). An isothermal column
 *  (20 C), no-flow side walls and base, fluid pressure fixed at 1 atm on the top boundary, air in the
 *  hydraulic-mixture mode (CVFEM_PHX_Scheme::ActivateAir(0., false, true, 0.)). The pore space holds as much air as
 *  the scheme allows: the pore water is at the floor of 2 % of the pore volume that NaClH2OPropertiesVisitorPHX
 *  always keeps for the H2O-NaCl fluid (min_h2o_nacl_fraction), and is immobile (below the residual liquid
 *  saturation 0.3). The air is the code's ideal gas, rho = p / (R_a T), so the exact solution is no flow and
 *      p(z) = p_top exp(-g (z - z_top) / (R_a T))   (z up),
 *  with R_a T = p_top / rho(T, p_top) taken from the code's own air equation of state (Air::Density) and g as in the
 *  scheme (9.80665 m s-2). The column is 2000 m tall, so the air density changes by 21 % from top to base.
 *
 *  Set-up: the exact profile is imposed as the initial pressure, the water-saturated state is computed from it
 *  (CVFEM_PHX_Scheme::InitialFluidPropertiesFromPTX), and CVFEM_AirPhaseLoader then replaces 98 % of the pore water
 *  by air of the exact density at every node. The top boundary keeps its pressure; without
 *  CVFEM_PHX_Scheme::OpenBoundaries no mass crosses it, so the model is closed for air. The run lasts 0.05 years in
 *  50 time steps of 0.001 years (at most max_steps_ = 400): any spurious air flow builds up within the first time
 *  steps, and the pressure diffuses over a cell (50 m) in about 4 days, so 18 days cover several cells.
 *
 *  Pass criteria, after the configured duration:
 *  - largest relative deviation |p - p_exact| / p_exact below max_pressure_error_ = 1e-4. The total pressure
 *    variation is 21 %; a density treated as constant over each cell, or taken from one end of every facet,
 *    deviates by about (dz / H_a)^2 / 2 = 1.7e-5 per cell (dz = 50 m, scale height H_a = R_a T / g = 8.6 km),
 *    i.e. up to 7e-4 over the 40 cells of the column, while a second-order discretisation deviates by
 *    O((dz / H_a)^2) = 3e-5 in total;
 *  - the air moved between interior control volumes (the transport scheme's "fluid flux integral", as a mean Darcy
 *    speed, as in test-205) below max_velocity_fraction_ = 1e-3 times the air buoyancy velocity scale
 *    k g rho_a / mu_a;
 *  - air mass conserved to max_air_mass_change_ = 1e-9 (relative): no air leaves a model without an open boundary,
 *    so the total air mass (sum of "fluid mass air" times "pore volume", computed here, not by the scheme) must not
 *    change beyond round-off.
 *  The transport at the fixed-pressure top nodes is reported for information.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_2D_VVCase and tests/CVFEM/CVFEM_HydrostaticColumn_Test.
 *
 *  Reference: the exponential pressure profile of an isothermal ideal-gas column is the barometric formula of an
 *  isothermal atmosphere, Wallace J.M. & Hobbs P.V. (2006), Atmospheric Science: An Introductory Survey, 2nd ed.,
 *  Academic Press, section 3.2 (the hydrostatic equation, scale height H = R T / g).
 *
 *  @date 20.11.2026
 */

#include "Test.h"

#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"
#include "ModelTime.h"
#include "CVFEM_AirPhaseLoader.h"

#include <array>
#include <memory>
#include <string>

namespace csmp {

/// Isothermal hydrostatic air column: the CVFEM PHX scheme must keep the exact exponential profile at rest.
class CVFEM_GasHydrostaticColumn_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt);
    /// facet_mean_gravity_density switches CVFEM_PHX_Scheme::UseFacetMeanGravityDensity.
    CVFEM_GasHydrostaticColumn_Test(const std::string& mesh_name, const std::string& config_name,
                                    bool facet_mean_gravity_density);
    ~CVFEM_GasHydrostaticColumn_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, builds the CVFEM PHX scheme with air, permeability tensor and pore volumes.
    void CreateModel();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Top elevation, top pressure, temperature and R_a T from the configured state.
    void ReadTopState();
    /// Exact isothermal hydrostatic pressure at elevation y.
    double ExactPressure(double y) const;
    /// Exact pressure on every node, then the water-saturated fluid properties from it.
    void InitializeWaterSaturatedColumn();
    /// Replaces all pore water above the 2 % floor by air of the exact density.
    void LoadAirToTheWaterFloor();
    /// Largest time step, transient preparation; the first time step is the largest.
    void StartTransient(double largest_time_step);
    /// Steps to max_time, at most max_steps_ steps (counts steps_).
    void AdvanceToEndTime(double max_time);
    /// Largest |p - p_exact| / p_exact over the nodes.
    double MaxRelativePressureError() const;
    /// Largest time-averaged air transport speed (fluid flux integral / run_time as a Darcy speed through half the
    /// control-volume perimeter) over interior nodes [0] and over nodes with a fixed pressure [1].
    std::array<double, 2> MaxTransportSpeed(double run_time) const;
    /// Air buoyancy velocity scale k g rho_a,max / mu_a,min from the current nodal state.
    double AirVelocityScale() const;
    /// Smallest and largest nodal air saturation.
    std::array<double, 2> AirSaturationRange() const;
    /// Writes the pressure error, transport speeds, air mass change and their limits.
    void ReportResults(double pressure_error, double velocity_scale, const std::array<double, 2>& transport_speed,
                       double air_mass_change) const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    bool facet_mean_gravity_density_;
    uint64_t steps_;
    double max_pressure_error_;
    double max_velocity_fraction_;
    double max_air_mass_change_;
    double y_top_;
    double pressure_top_;
    double temperature_;
    double gas_constant_times_temperature_;   ///< R_a T = p / rho_a [m2 s-2]

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    std::unique_ptr<CVFEM_AirPhaseLoader> air_phase_loader_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr uint64_t max_steps_{400U};   ///< bounds the run time if the time step collapses
    static constexpr double gravity_{9.80665};             ///< as CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV
    static constexpr double water_floor_saturation_{0.02};   ///< NaClH2OPropertiesVisitorPHX min_h2o_nacl_fraction
};

} // namespace csmp

#endif // CVFEM_GAS_HYDROSTATIC_COLUMN_TEST_H
