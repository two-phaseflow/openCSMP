// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_AIR_COMPRESSION_BOX_TEST_H
#define CVFEM_AIR_COMPRESSION_BOX_TEST_H

/** @file CVFEM_AirCompressionBox_Test.h
 *  @author Samuel Scott
 *  @brief Air storage of the CVFEM PHX scheme: excess air in a closed, isothermal porous box must equilibrate to the
 *  ideal-gas pressure while conserving the air and water mass.
 *  @details Closed box (all boundaries no-flow and insulated), isothermal (20 C), no gravity, pore water plus air in
 *  the hydraulic-mixture mode (CVFEM_PHX_Scheme::ActivateAir(0., false, true, 0.)). Initial state: uniform pressure
 *  p0 and liquid saturation s_l, and the air that fills the rest of the pore space at p0, (1 - s_l) rho_a(T, p0) per
 *  unit pore volume, times 1 + e (1 + cos(pi xi)), xi = (x - x_min) / (x_max - x_min) (CVFEM_AirPhaseLoader). With
 *  e = 0.25 the left wall holds 1.5 times the filling amount, the right wall the filling amount, the box 1.25 times
 *  on average. The excess
 *  air is the only disturbance, a smooth one (the first mode of the box, so the time steps are not limited by a sharp
 *  front): the left side is overpressured, air (and, above the residual liquid saturation 0.3, water) flows to the
 *  right until the pressure is uniform. Without capillary pressure and gravity any saturation distribution at uniform
 *  pressure is at rest, so the end state is fixed by the conserved masses alone: with the code's ideal-gas air,
 *  rho_a = p / (R_a T),
 *      p_end (V_pore - M_w / rho_l) = M_a R_a T,
 *  with M_a, M_w the air and water mass of the box, V_pore its pore volume and rho_l the liquid density at p_end
 *  (from the run). R_a T = p0 / rho_a(T, p0) comes from the code's own air equation of state (Air::Density).
 *
 *  Cases, each for one regime of the air storage term of the pressure equation (the cross term of
 *  NaClH2OPropertiesVisitorPHX::ComputeTotalCompressibility):
 *  - test-209a, s_l = 0.02, e = 0.25: pore water at the floor of 2 % of the pore volume that
 *    NaClH2OPropertiesVisitorPHX always keeps for the H2O-NaCl fluid (min_h2o_nacl_fraction); the excess air does
 *    not fit at p0 and is held back by the cap on the air saturation, 0.98 (the "clamped" regime, storage
 *    rho_a s_a phi / p);
 *  - test-209b, s_l = 0.5, e = 0.02: substantial, mobile water that backfills the pore volume the air gives up
 *    (storage rho_l s_a phi / p); the small excess keeps the time step, limited by the air velocity, large against
 *    the slow pressure diffusion;
 *  In 209a the water is immobile (below the residual liquid saturation 0.3).
 *
 *  Pass criteria, at the end of the configured duration:
 *  - total air mass and total water mass (sums of nodal mass times pore volume, computed here) conserved to
 *    max_mass_change_ = 1e-6 (relative): the box is closed;
 *  - pressure uniform: (p_max - p_min) / p_mean below max_pressure_spread_ = 1e-3;
 *  - p_mean within max_pressure_error_ = 1e-3 (relative) of p_end above, and p_mean V_a within the same of
 *    M_a R_a T, V_a the air volume of the scheme (sum of "saturation air" times pore volume): air held back by the
 *    saturation cap would show as V_a < M_a R_a T / p;
 *  - at every node and after every time step: s_l + s_v + s_a + s_h = 1 within max_closure_error_ = 1e-9, all
 *    saturations, the air mass and the pressure non-negative;
 *  - the run reaches the configured duration within max_steps_ = 400 time steps, without a collapse of the time
 *    step. The durations are 5 to 11 decay times (L / pi)^2 / D of the cosine mode (L = 200 m): 0.2 years in 209a
 *    and 209c (air diffusivity D = k p / (mu_a phi s_a) = 7e-3 m2 s-1, decay time 0.018 years), 7.5 years in 209b
 *    (D = k (kr_l rho_l / mu_l + kr_a rho_a / mu_a) / (rho_l s_a phi / p) = 8.6e-5 m2 s-1, decay time 1.5 years).
 *  The tolerances of 1e-3 leave room for a residual pressure gradient and for one known approximation: at the
 *  2 % floor the air saturation is capped at 0.98, so when the water is compressed the air cannot take up the
 *  volume it frees; about (0.02 / 0.98) (rho_l / rho_a) beta_l (p - p0) of the air, 1.5e-4 at p - p0 = 0.25 bar, is
 *  then carried as excess mass above the cap, and p V_a falls short of M_a R_a T by as much.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_2D_VVCase and tests/CVFEM/CVFEM_HydrostaticColumn_Test.
 *
 *  Reference: ideal-gas law (Boyle's law at constant temperature), e.g. Wallace J.M. & Hobbs P.V. (2006),
 *  Atmospheric Science: An Introductory Survey, 2nd ed., Academic Press, section 3.1.
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

#include <memory>
#include <string>

namespace csmp {

/// Closed, isothermal porous box with excess air: equilibrium pressure versus the ideal gas, mass conservation.
class CVFEM_AirCompressionBox_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt); liquid_saturation
    /// s_l of the whole box; excess_air e: mean relative excess of the air mass, see LoadAirWithSmoothExcess.
    CVFEM_AirCompressionBox_Test(const std::string& mesh_name, const std::string& config_name,
                                 double liquid_saturation, double excess_air);
    ~CVFEM_AirCompressionBox_Test() override;

    void run() override;

  private:
    /// Totals over the nodes, weighted by pore volume.
    struct BoxState {
      double air_mass{0.};         ///< [kg]
      double water_mass{0.};       ///< H2O-NaCl mass [kg]
      double pore_volume{0.};      ///< [m3]
      double air_volume{0.};       ///< sum of saturation air times pore volume [m3]
      double mean_pressure{0.};    ///< pore-volume weighted [Pa]
      double pressure_spread{0.};  ///< (p_max - p_min) / mean_pressure
      double liquid_density{0.};   ///< pore-volume weighted [kg m-3]
      double temperature{0.};      ///< pore-volume weighted [C]
    };

    /// Reads mesh and configuration, builds the CVFEM PHX scheme with air, permeability tensor and pore volumes.
    void CreateModel();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Water-saturated state at the configured uniform pressure p0 and temperature.
    void InitializeWaterSaturatedBox();
    /// Liquid saturation s_l everywhere; the air filling the rest of the pore space at p0, times
    /// 1 + e (1 + cos(pi xi)), xi = (x - x_min) / (x_max - x_min).
    void LoadAirWithSmoothExcess();
    /// Largest time step, transient preparation; the first time step is the largest.
    void StartTransient(double largest_time_step);
    /// Steps to max_time, at most max_steps_ steps, checking saturation closure and positivity after every step
    /// (counts steps_).
    void AdvanceToEndTime(double max_time);
    /// Worst saturation-closure error and smallest saturation, air mass or pressure over the nodes, accumulated in
    /// max_closure_error_found_ and min_value_found_.
    void CheckClosureAndPositivity();
    /// Masses, volumes, mean pressure and its spread, liquid density and temperature of the box.
    BoxState ComputeBoxState() const;
    /// R_a T = p0 / rho_a(T, p0) from the code's air equation of state [m2 s-2].
    double GasConstantTimesTemperature(double temperature) const;
    /// p_end from the conserved masses: M_a R_a T / (V_pore - M_w / rho_l).
    double IdealGasEquilibriumPressure(const BoxState& initial, const BoxState& final) const;
    /// p_mean V_a / (M_a R_a T) - 1 with the air volume V_a of the scheme.
    double IdealGasLawError(const BoxState& final) const;
    /// Writes the measured quantities and their limits.
    void ReportResults(const BoxState& initial, const BoxState& final) const;
    /// Pass/fail checks.
    void CheckResults(const BoxState& initial, const BoxState& final, double max_time);

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    double liquid_saturation_;
    double excess_air_;
    uint64_t steps_;
    double max_mass_change_;
    double max_pressure_spread_;
    double max_pressure_error_;
    double max_closure_error_;
    double initial_pressure_;
    double max_closure_error_found_;
    double min_value_found_;                  ///< smallest saturation, air mass or pressure during the run

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    std::unique_ptr<CVFEM_AirPhaseLoader> air_phase_loader_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr uint64_t max_steps_{400U};   ///< bounds the run time if the time step collapses
};

} // namespace csmp

#endif // CVFEM_AIR_COMPRESSION_BOX_TEST_H
