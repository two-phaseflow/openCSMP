// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_AirCompressionBox_Test.cpp
 *  @author Samuel Scott
 *  @brief Air storage of the CVFEM PHX scheme in a closed box, see CVFEM_AirCompressionBox_Test.h.
 *  @date 20.11.2026
 */

#include "CVFEM_AirCompressionBox_Test.h"

#include "Air.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace csmp {

CVFEM_AirCompressionBox_Test::CVFEM_AirCompressionBox_Test(const std::string& mesh_name,
                                                           const std::string& config_name,
                                                           double liquid_saturation, double excess_air)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    liquid_saturation_(liquid_saturation),
    excess_air_(excess_air),
    steps_(0U),
    max_mass_change_(1.e-6),
    max_pressure_spread_(1.e-3),
    max_pressure_error_(1.e-3),
    max_closure_error_(1.e-9),
    initial_pressure_(0.),
    max_closure_error_found_(0.),
    min_value_found_(std::numeric_limits<double>::max()),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    air_phase_loader_(nullptr),
    run_settings_()
{
}


CVFEM_AirCompressionBox_Test::~CVFEM_AirCompressionBox_Test() = default;


void CVFEM_AirCompressionBox_Test::CreateModel()
{
  ModelBuilder<2U> builder;
  builder.CreateVariablesFile(false, false, false, false, false, false);
  builder.CreateLimiterLogFile();

  model_ = std::make_unique<ANSYS_Model2D>(mesh_name_.c_str(), mesh_name_.c_str(), vars_name_.c_str());
  model_->InstantiateFiniteVolumes();
  property_database_ = &(model_->Database());

  InputDataManager<2U> model_configuration;
  model_configuration.ConfigureFromFile(*model_, config_name_.c_str(),
                                        false,   ///< region name from parameter range
                                        true,    ///< default property values
                                        true,    ///< regional property values
                                        true,    ///< boundary conditions for box-shaped model
                                        true,    ///< regional property conditions
                                        true,    ///< boundary conditions for arbitrary-shaped model
                                        run_settings_);

  const bool with_gravity{static_cast<bool>(model_->Read(property_database_->StorageKey("gravity flag")))};
  const std::vector<std::string> wells_list{};
  const std::map<std::string, WellConfiguration> well_configurations{};
  cvfem_phx_ = std::make_unique<CVFEM_PHX_Scheme<2U>>(*model_, with_gravity, true, wells_list, false, false,
                                                      SolverKind::PETSc, SolverKind::PETSc, well_configurations);
  const double heat_capacity_rock{model_->Region("Model").NodeVector().front()->Read(
      property_database_->StorageKey("nodal heat capacity rock"))};
  cvfem_phx_->ActivateAir(0., false, true, 0.);   // dry air, hydraulic gas mixture, no rain
  cvfem_phx_->SetRockHeatCapacity(heat_capacity_rock);
  cvfem_phx_->AddFluidContributionToHeatCapacity(true);

  permeability_visitor_ = std::make_unique<PermeabilityVisitor<2U>>(*model_);
  permeability_visitor_->AssignPermeabilityTensor(*model_);

  PoreVolumeVisitor<2U> pore_volume_visitor(*model_, "nodal porosity", "bulk volume", "pore volume", "thickness");
  model_->Accept(pore_volume_visitor);

  air_phase_loader_ = std::make_unique<CVFEM_AirPhaseLoader>(*model_, heat_capacity_rock);
}


void CVFEM_AirCompressionBox_Test::InitializeAllNaNsToZero()
{
  double minimum{0.};
  double maximum{0.};
  for (auto it = model_->Database().Begin(); it != model_->Database().End(); ++it) {
    model_->MinMaxOf((*it).first.c_str(), minimum, maximum);
    if (!std::isnan(minimum) && !std::isnan(maximum))
      continue;
    if ((*it).second.key.type == SCALAR)
      model_->InputPropertyValue((*it).first.c_str(), ScalarVariable(ANY, 0.));
    else if ((*it).second.key.type == VECTOR)
      model_->InputPropertyValue((*it).first.c_str(), VectorVariable<2U>(ANY, 0.));
    else if ((*it).second.key.type == TENSOR)
      model_->InputPropertyValue((*it).first.c_str(), TensorVariable<2U>(ANY, 0.));
  }
}


void CVFEM_AirCompressionBox_Test::InitializeWaterSaturatedBox()
{
  initial_pressure_ = model_->Region("Model").NodeVector().front()->Read(
      property_database_->StorageKey("fluid pressure"));
  cvfem_phx_->InitialFluidPropertiesFromPTX();
}


void CVFEM_AirCompressionBox_Test::LoadAirWithSmoothExcess()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  double x_min{std::numeric_limits<double>::max()};
  double x_max{std::numeric_limits<double>::lowest()};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    x_min = std::min(x_min, (*node_it)->x());
    x_max = std::max(x_max, (*node_it)->x());
  }
  const double pi{std::acos(-1.)};
  Air air;
  air_phase_loader_->Load(liquid_saturation_, [&](const Node<2U>& node) {
    const double filling_air{(1. - liquid_saturation_) * air.Density(node.Read(temperature_key), initial_pressure_)};
    const double excess{excess_air_ * (1. + std::cos(pi * (node.x() - x_min) / (x_max - x_min)))};
    return (1. + excess) * filling_air;
  });
}


void CVFEM_AirCompressionBox_Test::StartTransient(double largest_time_step)
{
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(largest_time_step);
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(largest_time_step);
}


void CVFEM_AirCompressionBox_Test::CheckClosureAndPositivity()
{
  const std::vector<Index> saturation_keys{property_database_->StorageKey("saturation liquid"),
                                           property_database_->StorageKey("saturation vapor"),
                                           property_database_->StorageKey("saturation air"),
                                           property_database_->StorageKey("saturation halite")};
  const Index air_mass_key(property_database_->StorageKey("fluid mass air"));
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    double saturation_sum{0.};
    for (const Index& key : saturation_keys) {
      saturation_sum += (*node_it)->Read(key);
      min_value_found_ = std::min(min_value_found_, (*node_it)->Read(key));
    }
    max_closure_error_found_ = std::max(max_closure_error_found_, std::fabs(saturation_sum - 1.));
    min_value_found_ = std::min({min_value_found_, (*node_it)->Read(air_mass_key), (*node_it)->Read(pressure_key)});
  }
}


void CVFEM_AirCompressionBox_Test::AdvanceToEndTime(double max_time)
{
  double& model_time(ModelTime::Instance().modelTime);
  steps_ = 0U;
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps_;
    CheckClosureAndPositivity();
    if (steps_ <= 5U || steps_ % 50U == 0U) {
      const BoxState state{ComputeBoxState()};
      std::cerr << "CVFEM_AirCompressionBox_Test::AdvanceToEndTime: step " << steps_ << ", " << model_time / year_
                << " yr, mean pressure " << state.mean_pressure << " Pa, spread " << state.pressure_spread
                << std::endl;
    }
    if (dt < 1. || steps_ >= max_steps_)
      break;
  }
}


CVFEM_AirCompressionBox_Test::BoxState CVFEM_AirCompressionBox_Test::ComputeBoxState() const
{
  const Index air_mass_key(property_database_->StorageKey("fluid mass air"));
  const Index total_mass_key(property_database_->StorageKey("fluid density"));
  const Index pore_volume_key(property_database_->StorageKey("pore volume"));
  const Index air_saturation_key(property_database_->StorageKey("saturation air"));
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index liquid_density_key(property_database_->StorageKey("density liquid"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  BoxState state;
  double pressure_min{std::numeric_limits<double>::max()};
  double pressure_max{std::numeric_limits<double>::lowest()};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const Node<2U>& node(**node_it);
    const double pore_volume{node.Read(pore_volume_key)};
    const double pressure{node.Read(pressure_key)};
    state.air_mass += node.Read(air_mass_key) * pore_volume;
    state.water_mass += (node.Read(total_mass_key) - node.Read(air_mass_key)) * pore_volume;
    state.pore_volume += pore_volume;
    state.air_volume += node.Read(air_saturation_key) * pore_volume;
    state.mean_pressure += pressure * pore_volume;
    state.liquid_density += node.Read(liquid_density_key) * pore_volume;
    state.temperature += node.Read(temperature_key) * pore_volume;
    pressure_min = std::min(pressure_min, pressure);
    pressure_max = std::max(pressure_max, pressure);
  }
  state.mean_pressure /= state.pore_volume;
  state.liquid_density /= state.pore_volume;
  state.temperature /= state.pore_volume;
  state.pressure_spread = (pressure_max - pressure_min) / state.mean_pressure;
  return state;
}


double CVFEM_AirCompressionBox_Test::GasConstantTimesTemperature(double temperature) const
{
  Air air;
  return initial_pressure_ / air.Density(temperature, initial_pressure_);
}


double CVFEM_AirCompressionBox_Test::IdealGasEquilibriumPressure(const BoxState& initial, const BoxState& final) const
{
  return initial.air_mass * GasConstantTimesTemperature(final.temperature)
         / (initial.pore_volume - initial.water_mass / final.liquid_density);
}


double CVFEM_AirCompressionBox_Test::IdealGasLawError(const BoxState& final) const
{
  return final.mean_pressure * final.air_volume / (final.air_mass * GasConstantTimesTemperature(final.temperature))
         - 1.;
}


void CVFEM_AirCompressionBox_Test::ReportResults(const BoxState& initial, const BoxState& final) const
{
  const double equilibrium_pressure{IdealGasEquilibriumPressure(initial, final)};
  std::cerr << "CVFEM_AirCompressionBox_Test::run: " << mesh_name_ << ", liquid saturation " << liquid_saturation_
            << ", excess air " << excess_air_ << ": " << steps_ << " steps to "
            << ModelTime::Instance().modelTime / year_ << " yr" << std::endl
            << "  relative change of the total air mass "
            << (final.air_mass - initial.air_mass) / initial.air_mass << ", water mass "
            << (final.water_mass - initial.water_mass) / initial.water_mass << " (limit " << max_mass_change_ << ")"
            << std::endl
            << "  pressure: initial " << initial_pressure_ << " Pa, mean " << final.mean_pressure
            << " Pa, ideal-gas equilibrium " << equilibrium_pressure << " Pa (relative error "
            << final.mean_pressure / equilibrium_pressure - 1. << ", limit " << max_pressure_error_ << "), spread "
            << final.pressure_spread << " (limit " << max_pressure_spread_ << ")" << std::endl
            << "  p V_a / (M_a R_a T) - 1 = " << IdealGasLawError(final) << " (limit " << max_pressure_error_
            << "), air volume fraction of the pore space " << final.air_volume / final.pore_volume << std::endl
            << "  worst saturation closure error " << max_closure_error_found_ << " (limit " << max_closure_error_
            << "), smallest saturation, air mass or pressure " << min_value_found_ << std::endl;
}


void CVFEM_AirCompressionBox_Test::CheckResults(const BoxState& initial, const BoxState& final, double max_time)
{
  const double equilibrium_pressure{IdealGasEquilibriumPressure(initial, final)};
  do_test(ModelTime::Instance().modelTime >= max_time, "Air compression box: ran the configured duration", __FILE__,
          __LINE__);
  do_test(std::fabs(final.air_mass - initial.air_mass) < max_mass_change_ * initial.air_mass,
          "Air compression box: total air mass conserved", __FILE__, __LINE__);
  do_test(std::fabs(final.water_mass - initial.water_mass) < max_mass_change_ * initial.water_mass,
          "Air compression box: total water mass conserved", __FILE__, __LINE__);
  do_test(final.pressure_spread < max_pressure_spread_, "Air compression box: pressure equilibrated", __FILE__,
          __LINE__);
  do_test(std::fabs(final.mean_pressure / equilibrium_pressure - 1.) < max_pressure_error_,
          "Air compression box: pressure equals the ideal-gas equilibrium pressure of the conserved masses", __FILE__,
          __LINE__);
  do_test(std::fabs(IdealGasLawError(final)) < max_pressure_error_,
          "Air compression box: p V_a = M_a R_a T with the air volume of the scheme", __FILE__, __LINE__);
  do_test(max_closure_error_found_ < max_closure_error_, "Air compression box: saturations sum to one", __FILE__,
          __LINE__);
  do_test(min_value_found_ >= 0., "Air compression box: saturations, air mass and pressure non-negative", __FILE__,
          __LINE__);
}


void CVFEM_AirCompressionBox_Test::run()
{
  CreateModel();
  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};

  InitializeAllNaNsToZero();
  InitializeWaterSaturatedBox();
  LoadAirWithSmoothExcess();
  CheckClosureAndPositivity();
  const BoxState initial{ComputeBoxState()};
  StartTransient(run_settings_.TimeIncrement() * year_);
  AdvanceToEndTime(max_time);

  const BoxState final{ComputeBoxState()};
  ReportResults(initial, final);
  CheckResults(initial, final, max_time);
}

} // namespace csmp
