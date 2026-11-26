// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_RechargeMound_Test.cpp
 *  @author Samuel Scott
 *  @brief Steady groundwater mound under uniform rain against the Dupuit solution, see CVFEM_RechargeMound_Test.h.
 *  @date 20.11.2026
 */

#include "CVFEM_RechargeMound_Test.h"

#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "Air.h"
#include "PoreVolumeVisitor.h"
#include "WellConfigurationFile.h"
#include "ModelTime.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

namespace csmp {

CVFEM_RechargeMound_Test::CVFEM_RechargeMound_Test(const std::string& mesh_name, const std::string& config_name,
                                                   double rain_mm_per_year, double outlet_water_level)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    rain_mm_per_year_(rain_mm_per_year),
    outlet_water_level_(outlet_water_level),
    steps_(0U),
    x_min_(0.),
    x_max_(0.),
    y_base_(0.),
    y_top_(0.),
    max_divide_error_(0.02),
    max_column_error_(0.05),
    vadose_saturation_(0.3),
    air_weight_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    air_phase_loader_(nullptr),
    run_settings_()
{
}


CVFEM_RechargeMound_Test::~CVFEM_RechargeMound_Test() = default;


void CVFEM_RechargeMound_Test::CreateModel()
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
  cvfem_phx_->SetRockHeatCapacity(heat_capacity_rock);
  cvfem_phx_->AddFluidContributionToHeatCapacity(true);
  cvfem_phx_->ActivateAir(0., false, true, rain_mm_per_year_);   // dry air, hydraulic gas mixture, rain
  cvfem_phx_->OpenBoundaries(0.);                                 // fresh water flows in and out at fixed-T nodes

  permeability_visitor_ = std::make_unique<PermeabilityVisitor<2U>>(*model_);
  permeability_visitor_->AssignPermeabilityTensor(*model_);

  PoreVolumeVisitor<2U> pore_volume_visitor(*model_, "nodal porosity", "bulk volume", "pore volume", "thickness");
  model_->Accept(pore_volume_visitor);
  air_phase_loader_ = std::make_unique<CVFEM_AirPhaseLoader>(*model_, heat_capacity_rock);

  x_min_ = y_base_ = std::numeric_limits<double>::max();
  x_max_ = y_top_ = std::numeric_limits<double>::lowest();
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    x_min_ = std::min(x_min_, (*node_it)->x());
    x_max_ = std::max(x_max_, (*node_it)->x());
    y_base_ = std::min(y_base_, (*node_it)->y());
    y_top_ = std::max(y_top_, (*node_it)->y());
  }
}


void CVFEM_RechargeMound_Test::InitializeAllNaNsToZero()
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


void CVFEM_RechargeMound_Test::InitializeFlatWaterTable()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const double water_table{y_base_ + outlet_water_level_};
  // air at rest above the water table (from p_atm at the land surface), water at rest below it
  Air air;
  const double air_density{air.Density(model_->Region("Model").NodeVector().front()->Read(temperature_key),
                                       atmospheric_pressure_)};
  air_weight_ = air_density * gravity_;
  const double water_table_pressure{atmospheric_pressure_ + air_density * gravity_ * (y_top_ - water_table)};
  for (uint32_t pass = 0U; pass < 2U; ++pass) {
    const Region<2U>& region(model_->Region("Model"));
    for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
      ScalarVariable pressure;
      (*node_it)->Read(pressure_key, pressure);
      const bool outlet{std::fabs((*node_it)->x() - x_max_) < 1.e-6};
      if (pressure.Flag() == DIRICH && !outlet)
        continue;   // land surface
      const double density{pass == 0U ? 1000. : (*node_it)->Read(density_key)};
      const double y{(*node_it)->y()};
      pressure() = y >= water_table ? atmospheric_pressure_ + air_density * gravity_ * (y_top_ - y)
                                    : water_table_pressure + density * gravity_ * (water_table - y);
      if (outlet)
        pressure.Flag() = DIRICH;
      (*node_it)->Store(pressure_key, pressure);
    }
    cvfem_phx_->InitialFluidPropertiesFromPTX();
  }
}


void CVFEM_RechargeMound_Test::LoadAirAboveWaterTable()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const auto saturation = [&](const Node<2U>& node) {
    return node.y() < y_base_ + outlet_water_level_ - 1.e-6 ? 1. : vadose_saturation_;
  };
  Air air;
  air_phase_loader_->Load(saturation, [&](const Node<2U>& node) {
    return (1. - saturation(node)) * air.Density(node.Read(temperature_key), node.Read(pressure_key));
  });
}


void CVFEM_RechargeMound_Test::Advance(double max_time)
{
  double& model_time(ModelTime::Instance().modelTime);
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps_;
    if (steps_ <= 10U || steps_ % 50U == 0U)
    {
      // smallest CFL-allowed time step per phase
      double courant_min_liquid{0.}, courant_max{0.}, courant_min_air{0.};
      model_->MinMaxOf("courant liquid", courant_min_liquid, courant_max);
      model_->MinMaxOf("courant air", courant_min_air, courant_max);
      std::cerr << "CVFEM_RechargeMound_Test::Advance: step " << steps_ << ", " << model_time / year_ << " yr, dt "
                << dt / year_ << " yr, water table at the divide " << WaterTable().begin()->second
                << " m, CFL dt liquid " << courant_min_liquid / year_ << " yr, air " << courant_min_air / year_ << " yr";
      std::cerr << std::endl;
    }
    if (dt < 1.)
      break;
  }
}


std::map<double, double> CVFEM_RechargeMound_Test::WaterTable() const
{
  // per column: the highest node of the saturated zone (no air) and the node below it; the water table is where their
  // pressure gradient meets the pressure of air at rest, p_atm + rho_air g (y_top - y)
  constexpr double saturated{1.e-3};   // largest air saturation of a saturated node
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index air_saturation_key(property_database_->StorageKey("saturation air"));
  const double air_weight{air_weight_};
  std::map<double, std::map<double, std::pair<double, double>>> columns;   // x -> (y -> (p, s_air))
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    columns[std::round((*node_it)->x() * 1.e3) / 1.e3][(*node_it)->y()] = {(*node_it)->Read(pressure_key),
                                                                          (*node_it)->Read(air_saturation_key)};
  std::map<double, double> water_table;
  for (const auto& [x, column] : columns) {
    if (column.size() < 3U)   // a node off the mesh columns (test-211c shifts one node row sideways), not a column
      continue;
    auto highest{column.begin()};
    while (std::next(highest) != column.end() && std::next(highest)->second.second < saturated)
      ++highest;
    if (highest == column.begin() || highest->second.second >= saturated) {
      water_table[x] = 0.;
      continue;
    }
    const auto below{std::prev(highest)};
    const double gradient{(highest->second.first - below->second.first) / (highest->first - below->first)};   // < 0
    // p_h + gradient (y - y_h) = p_atm + air_weight (y_top - y)
    const double y{(atmospheric_pressure_ + air_weight * y_top_ - highest->second.first + gradient * highest->first)
                   / (gradient + air_weight)};
    water_table[x] = y - y_base_;
  }
  return water_table;
}


double CVFEM_RechargeMound_Test::RainMassRate() const
{
  const Index rain_key(property_database_->StorageKey("rain recharge top"));
  const Index pore_volume_key(property_database_->StorageKey("pore volume"));
  double rate{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    rate += (*node_it)->Read(rain_key) * (*node_it)->Read(pore_volume_key);
  return rate;
}


void CVFEM_RechargeMound_Test::run()
{
  CreateModel();
  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};
  InitializeAllNaNsToZero();
  // Darcy-flux CFL criterion (UpwindControlVisitor, option B): the pore-velocity criterion divides by the phase
  // saturation, so the trace of air at the rising water table would set the time step
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), false);
  cvfem_phx_->SetLargestTimeStep(run_settings_.TimeIncrement() * year_);
  InitializeFlatWaterTable();
  LoadAirAboveWaterTable();
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(run_settings_.TimeIncrement() * year_ / 10.);

  Advance(0.9 * max_time);
  const double divide_before{WaterTable().begin()->second};
  Advance(max_time);
  const std::map<double, double> water_table{WaterTable()};

  // Dupuit: recharge from the assigned rain mass, conductivity from the liquid at the base of the divide
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Index viscosity_key(property_database_->StorageKey("viscosity liquid"));
  const Node<2U>* base_node{nullptr};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if (std::fabs((*node_it)->x() - x_min_) < 1.e-6 && std::fabs((*node_it)->y() - y_base_) < 1.e-6)
      base_node = *node_it;
  const double density{base_node->Read(density_key)};
  const double length{x_max_ - x_min_};
  const double recharge{RainMassRate() / (length * density)};   // [m/s] groundwater volume per unit area
  const double permeability{(*region.CellsBegin())->Read(property_database_->StorageKey("permeability"))};
  const double conductivity{permeability * density * gravity_ / base_node->Read(viscosity_key)};
  auto dupuit = [&](double x) {
    return std::sqrt(outlet_water_level_ * outlet_water_level_
                     + recharge / conductivity * (length * length - (x - x_min_) * (x - x_min_)));
  };

  std::ofstream csv(config_name_ + "-water-table.csv");
  csv << "x,water_table,dupuit" << std::endl;
  double max_column_deviation{0.};
  for (const auto& [x, h] : water_table) {
    csv << x << "," << h << "," << dupuit(x) << std::endl;
    if (x < x_max_ - 1.e-6 && std::next(water_table.find(x)) != water_table.end()
        && std::next(water_table.find(x))->first < x_max_ - 1.e-6)   // all but the last two columns (outlet)
      max_column_deviation = std::max(max_column_deviation, std::fabs(h - dupuit(x)) / dupuit(x));
  }
  const double divide{water_table.begin()->second};
  const double divide_deviation{std::fabs(divide - dupuit(x_min_)) / dupuit(x_min_)};

  std::cerr << "CVFEM_RechargeMound_Test::run: " << mesh_name_ << ": " << steps_ << " steps to "
            << model_time / year_ << " yr" << std::endl
            << "  recharge W " << recharge * year_ << " m/yr, conductivity K " << conductivity << " m/s, W/K "
            << recharge / conductivity << ", outlet water level " << outlet_water_level_ << " m above the base"
            << std::endl
            << "  water table at the divide " << divide << " m, Dupuit " << dupuit(x_min_) << " m (relative deviation "
            << divide_deviation << ", limit " << max_divide_error_ << "); change over the last 10 % of the run "
            << divide - divide_before << " m" << std::endl
            << "  largest relative deviation over the columns (outlet excluded) " << max_column_deviation
            << " (limit " << max_column_error_ << ")" << std::endl;

  do_test(model_time >= max_time, "Recharge mound: ran the configured duration", __FILE__, __LINE__);
  do_test(divide_deviation < max_divide_error_, "Recharge mound: water table at the divide matches Dupuit",
          __FILE__, __LINE__);
  do_test(max_column_deviation < max_column_error_, "Recharge mound: water table profile matches Dupuit",
          __FILE__, __LINE__);
}

} // namespace csmp
