// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_GasHydrostaticColumn_Test.cpp
 *  @author Samuel Scott
 *  @brief Well-balancedness of the air phase of the CVFEM PHX scheme, see CVFEM_GasHydrostaticColumn_Test.h.
 *  @date 20.11.2026
 */

#include "CVFEM_GasHydrostaticColumn_Test.h"

#include "Air.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace csmp {

CVFEM_GasHydrostaticColumn_Test::CVFEM_GasHydrostaticColumn_Test(const std::string& mesh_name,
                                                                 const std::string& config_name,
                                                                 bool facet_mean_gravity_density)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    facet_mean_gravity_density_(facet_mean_gravity_density),
    steps_(0U),
    max_pressure_error_(1.e-4),
    max_velocity_fraction_(1.e-3),
    max_air_mass_change_(1.e-9),
    y_top_(0.),
    pressure_top_(0.),
    temperature_(0.),
    gas_constant_times_temperature_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    air_phase_loader_(nullptr),
    run_settings_()
{
}


CVFEM_GasHydrostaticColumn_Test::~CVFEM_GasHydrostaticColumn_Test() = default;


void CVFEM_GasHydrostaticColumn_Test::CreateModel()
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
  cvfem_phx_->UseFacetMeanGravityDensity(facet_mean_gravity_density_);

  permeability_visitor_ = std::make_unique<PermeabilityVisitor<2U>>(*model_);
  permeability_visitor_->AssignPermeabilityTensor(*model_);

  PoreVolumeVisitor<2U> pore_volume_visitor(*model_, "nodal porosity", "bulk volume", "pore volume", "thickness");
  model_->Accept(pore_volume_visitor);

  air_phase_loader_ = std::make_unique<CVFEM_AirPhaseLoader>(*model_, heat_capacity_rock);
}


void CVFEM_GasHydrostaticColumn_Test::InitializeAllNaNsToZero()
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


void CVFEM_GasHydrostaticColumn_Test::ReadTopState()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  const Node<2U>* top_node{region.NodeVector().front()};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if ((*node_it)->y() > top_node->y())
      top_node = *node_it;
  y_top_ = top_node->y();
  pressure_top_ = top_node->Read(pressure_key);
  temperature_ = top_node->Read(temperature_key);
  Air air;
  gas_constant_times_temperature_ = pressure_top_ / air.Density(temperature_, pressure_top_);
}


double CVFEM_GasHydrostaticColumn_Test::ExactPressure(double y) const
{
  return pressure_top_ * std::exp(-gravity_ * (y - y_top_) / gas_constant_times_temperature_);
}


void CVFEM_GasHydrostaticColumn_Test::InitializeWaterSaturatedColumn()
{
  ReadTopState();
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    ScalarVariable pressure;
    (*node_it)->Read(pressure_key, pressure);
    pressure() = ExactPressure((*node_it)->y());
    (*node_it)->Store(pressure_key, pressure);
  }
  cvfem_phx_->InitialFluidPropertiesFromPTX();
}


void CVFEM_GasHydrostaticColumn_Test::LoadAirToTheWaterFloor()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  Air air;
  air_phase_loader_->Load(water_floor_saturation_, [&](const Node<2U>& node) {
    return (1. - water_floor_saturation_) * air.Density(node.Read(temperature_key), node.Read(pressure_key));
  });
}


void CVFEM_GasHydrostaticColumn_Test::StartTransient(double largest_time_step)
{
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(largest_time_step);
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(largest_time_step);
}


void CVFEM_GasHydrostaticColumn_Test::AdvanceToEndTime(double max_time)
{
  double& model_time(ModelTime::Instance().modelTime);
  steps_ = 0U;
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps_;
    if (steps_ <= 5U || steps_ % 50U == 0U)
      std::cerr << "CVFEM_GasHydrostaticColumn_Test::AdvanceToEndTime: step " << steps_ << ", "
                << model_time / year_ << " yr, max relative pressure error " << MaxRelativePressureError()
                << std::endl;
    if (dt < 1. || steps_ >= max_steps_)
      break;
  }
}


double CVFEM_GasHydrostaticColumn_Test::MaxRelativePressureError() const
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  double maximum{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double exact{ExactPressure((*node_it)->y())};
    maximum = std::max(maximum, std::fabs((*node_it)->Read(pressure_key) - exact) / exact);
  }
  return maximum;
}


std::array<double, 2> CVFEM_GasHydrostaticColumn_Test::MaxTransportSpeed(double run_time) const
{
  const Index flux_integral_key(property_database_->StorageKey("fluid flux integral"));
  const Index density_key(property_database_->StorageKey("density air"));
  const Index bulk_volume_key(property_database_->StorageKey("bulk volume"));
  const Index thickness_key(property_database_->StorageKey("thickness"));
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  std::array<double, 2> maximum{0., 0.};   // [interior nodes, nodes with a fixed pressure]
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double thickness{(*node_it)->Read(thickness_key)};
    const double area{(*node_it)->Read(bulk_volume_key) / thickness};
    // mass that left the control volume (the liquid is immobile), as a mean Darcy speed through half its
    // perimeter (~ 2 sqrt(area))
    const double speed{(*node_it)->Read(flux_integral_key) / run_time
                       / ((*node_it)->Read(density_key) * 2. * std::sqrt(area) * thickness)};
    ScalarVariable pressure;
    (*node_it)->Read(pressure_key, pressure);
    const bool fixed_pressure{pressure.Flag() == DIRICH};
    maximum[fixed_pressure ? 1U : 0U] = std::max(maximum[fixed_pressure ? 1U : 0U], speed);
  }
  return maximum;
}


double CVFEM_GasHydrostaticColumn_Test::AirVelocityScale() const
{
  const Index density_key(property_database_->StorageKey("density air"));
  const Index viscosity_key(property_database_->StorageKey("viscosity air"));
  const Index permeability_key(property_database_->StorageKey("permeability"));
  const Region<2U>& region(model_->Region("Model"));
  double density_max{0.};
  double viscosity_min{std::numeric_limits<double>::max()};
  double permeability_max{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    density_max = std::max(density_max, (*node_it)->Read(density_key));
    viscosity_min = std::min(viscosity_min, (*node_it)->Read(viscosity_key));
  }
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it)
    permeability_max = std::max(permeability_max, (*element_it)->Read(permeability_key));
  return permeability_max * gravity_ * density_max / viscosity_min;
}


std::array<double, 2> CVFEM_GasHydrostaticColumn_Test::AirSaturationRange() const
{
  const Index saturation_key(property_database_->StorageKey("saturation air"));
  std::array<double, 2> range{std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    range[0] = std::min(range[0], (*node_it)->Read(saturation_key));
    range[1] = std::max(range[1], (*node_it)->Read(saturation_key));
  }
  return range;
}


void CVFEM_GasHydrostaticColumn_Test::ReportResults(double pressure_error, double velocity_scale,
                                                    const std::array<double, 2>& transport_speed,
                                                    double air_mass_change) const
{
  const std::array<double, 2> air_saturation{AirSaturationRange()};
  std::cerr << "CVFEM_GasHydrostaticColumn_Test::run: " << mesh_name_ << " (facet-mean gravity density "
            << (facet_mean_gravity_density_ ? "on" : "off") << "): " << steps_ << " steps to "
            << ModelTime::Instance().modelTime / year_ << " yr, R_a T " << gas_constant_times_temperature_
            << " m2 s-2, air velocity scale " << velocity_scale << " m/s" << std::endl
            << "  max relative pressure error " << pressure_error << " (limit " << max_pressure_error_ << ")"
            << std::endl
            << "  transported air as mean Darcy speed: interior nodes " << transport_speed[0] << " m/s (ratio "
            << transport_speed[0] / velocity_scale << ", limit " << max_velocity_fraction_ << "), fixed-pressure nodes "
            << transport_speed[1] << " m/s (ratio " << transport_speed[1] / velocity_scale << ")" << std::endl
            << "  relative change of the total air mass " << air_mass_change << " (limit " << max_air_mass_change_
            << "), air saturation " << air_saturation[0] << " to " << air_saturation[1] << std::endl;
}


void CVFEM_GasHydrostaticColumn_Test::run()
{
  CreateModel();
  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};

  InitializeAllNaNsToZero();
  InitializeWaterSaturatedColumn();
  LoadAirToTheWaterFloor();
  const double initial_air_mass{air_phase_loader_->TotalAirMass()};
  StartTransient(run_settings_.TimeIncrement() * year_);
  AdvanceToEndTime(max_time);

  const double pressure_error{MaxRelativePressureError()};
  const double velocity_scale{AirVelocityScale()};
  const std::array<double, 2> transport_speed{MaxTransportSpeed(model_time)};
  const double air_mass_change{std::fabs(air_phase_loader_->TotalAirMass() - initial_air_mass) / initial_air_mass};
  ReportResults(pressure_error, velocity_scale, transport_speed, air_mass_change);

  do_test(model_time >= max_time, "Gas hydrostatic column: ran the configured duration", __FILE__, __LINE__);
  do_test(pressure_error < max_pressure_error_,
          "Gas hydrostatic column: pressure within the limit of the exact isothermal profile", __FILE__, __LINE__);
  do_test(transport_speed[0] < max_velocity_fraction_ * velocity_scale,
          "Gas hydrostatic column: air transported between interior control volumes below the limit fraction of the "
          "air velocity scale", __FILE__, __LINE__);
  do_test(air_mass_change < max_air_mass_change_, "Gas hydrostatic column: total air mass conserved", __FILE__,
          __LINE__);
}

} // namespace csmp
