// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_AxisymmetricAnnulus_Test.cpp
 *  @author Samuel Scott
 *  @brief Steady radial heat transfer in an axisymmetric annulus, see CVFEM_AxisymmetricAnnulus_Test.h.
 *  @date 20.11.2026
 */

#include "CVFEM_AxisymmetricAnnulus_Test.h"

#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "PoreVolumeVisitor.h"
#include "WellConfigurationFile.h"
#include "ModelTime.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace csmp {

CVFEM_AxisymmetricAnnulus_Test::CVFEM_AxisymmetricAnnulus_Test(const std::string& mesh_name, const std::string& config_name)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    steps_(0U),
    inner_radius_(0.),
    outer_radius_(0.),
    height_(0.),
    inner_temperature_(0.),
    outer_temperature_(0.),
    peclet_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    run_settings_(),
    axisymmetric_{true, 0., 2. * M_PI}
{
}


CVFEM_AxisymmetricAnnulus_Test::~CVFEM_AxisymmetricAnnulus_Test() = default;


double CVFEM_AxisymmetricAnnulus_Test::CreateModel()
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
  // before the scheme and the pore volumes read it
  const double revolved_volume{setAxisymmetricThickness(*model_, axisymmetric_)};

  const bool with_gravity{static_cast<bool>(model_->Read(property_database_->StorageKey("gravity flag")))};
  const std::vector<std::string> wells_list{};
  const std::map<std::string, WellConfiguration> well_configurations{};
  cvfem_phx_ = std::make_unique<CVFEM_PHX_Scheme<2U>>(*model_, with_gravity, true, wells_list, false, false,
                                                      SolverKind::PETSc, SolverKind::PETSc, well_configurations);
  const double heat_capacity_rock{model_->Region("Model").NodeVector().front()->Read(
      property_database_->StorageKey("nodal heat capacity rock"))};
  cvfem_phx_->SetRockHeatCapacity(heat_capacity_rock);
  cvfem_phx_->AddFluidContributionToHeatCapacity(true);

  permeability_visitor_ = std::make_unique<PermeabilityVisitor<2U>>(*model_);
  permeability_visitor_->AssignPermeabilityTensor(*model_);

  PoreVolumeVisitor<2U> pore_volume_visitor(*model_, "nodal porosity", "bulk volume", "pore volume", "thickness");
  model_->Accept(pore_volume_visitor);
  return revolved_volume;
}


void CVFEM_AxisymmetricAnnulus_Test::InitializeAllNaNsToZero()
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


void CVFEM_AxisymmetricAnnulus_Test::ReadGeometryAndWallTemperatures()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  inner_radius_ = std::numeric_limits<double>::max();
  outer_radius_ = std::numeric_limits<double>::lowest();
  double y_bottom{std::numeric_limits<double>::max()};
  double y_top{std::numeric_limits<double>::lowest()};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    inner_radius_ = std::min(inner_radius_, (*node_it)->x());
    outer_radius_ = std::max(outer_radius_, (*node_it)->x());
    y_bottom = std::min(y_bottom, (*node_it)->y());
    y_top = std::max(y_top, (*node_it)->y());
  }
  height_ = y_top - y_bottom;
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    if (std::fabs((*node_it)->x() - inner_radius_) < 1.e-6) inner_temperature_ = (*node_it)->Read(temperature_key);
    if (std::fabs((*node_it)->x() - outer_radius_) < 1.e-6) outer_temperature_ = (*node_it)->Read(temperature_key);
  }
}


void CVFEM_AxisymmetricAnnulus_Test::Advance(double max_time)
{
  double& model_time(ModelTime::Instance().modelTime);
  steps_ = 0U;
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps_;
    if (dt < 1.)
      break;
  }
}


std::array<double, 2> CVFEM_AxisymmetricAnnulus_Test::RadialMassRate() const
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Index viscosity_key(property_database_->StorageKey("viscosity liquid"));
  const Index permeability_key(property_database_->StorageKey("permeability"));
  const Region<2U>& region(model_->Region("Model"));
  std::vector<double> mass_flux_times_radius;   // r rho q_r per element
  std::vector<double> areas;
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it) {
    const Element<2U>& element(**element_it);
    std::array<double, 3> x{};
    std::array<double, 3> y{};
    std::array<double, 3> p{};
    double density{0.};
    double viscosity{0.};
    for (uint32_t i = 0U; i < 3U; ++i) {
      x[i] = element.N(i)->x();
      y[i] = element.N(i)->y();
      p[i] = element.N(i)->Read(pressure_key);
      density += element.N(i)->Read(density_key) / 3.;
      viscosity += element.N(i)->Read(viscosity_key) / 3.;
    }
    // x-derivative of the linear pressure interpolant on the triangle
    const double twice_area{(x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])};
    const double dp_dx{((p[1] - p[0]) * (y[2] - y[0]) - (p[2] - p[0]) * (y[1] - y[0])) / twice_area};
    const double darcy_flux{-element.Read(permeability_key) / viscosity * dp_dx};
    mass_flux_times_radius.push_back(element.BaryCenter()[0] * density * darcy_flux);
    areas.push_back(0.5 * std::fabs(twice_area));
  }
  double sum{0.};
  double area{0.};
  for (size_t e = 0U; e < areas.size(); ++e) {
    sum += areas[e] * mass_flux_times_radius[e];
    area += areas[e];
  }
  const double mean{sum / area};
  double deviation{0.};
  for (double value : mass_flux_times_radius)
    deviation = std::max(deviation, std::fabs(value - mean));
  return {2. * M_PI * height_ * mean, mean != 0. ? deviation / std::fabs(mean) : 0.};
}


double CVFEM_AxisymmetricAnnulus_Test::WallEnthalpy(double x) const
{
  const Index enthalpy_key(property_database_->StorageKey("enthalpy liquid"));
  const Region<2U>& region(model_->Region("Model"));
  double sum{0.};
  uint32_t count{0U};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if (std::fabs((*node_it)->x() - x) < 1.e-6) {
      sum += (*node_it)->Read(enthalpy_key);
      ++count;
    }
  return sum / count;
}


double CVFEM_AxisymmetricAnnulus_Test::ExactTemperature(double r) const
{
  if (std::fabs(peclet_) < 1.e-6)
    return ConductionTemperature(r);
  return inner_temperature_ + (outer_temperature_ - inner_temperature_)
                                  * (std::pow(r, peclet_) - std::pow(inner_radius_, peclet_))
                                  / (std::pow(outer_radius_, peclet_) - std::pow(inner_radius_, peclet_));
}


double CVFEM_AxisymmetricAnnulus_Test::ConductionTemperature(double r) const
{
  return inner_temperature_ + (outer_temperature_ - inner_temperature_) * std::log(r / inner_radius_)
                                  / std::log(outer_radius_ / inner_radius_);
}


double CVFEM_AxisymmetricAnnulus_Test::PlaneTemperature(double x) const
{
  return inner_temperature_ + (outer_temperature_ - inner_temperature_) * (x - inner_radius_)
                                  / (outer_radius_ - inner_radius_);
}


double CVFEM_AxisymmetricAnnulus_Test::MaxTemperatureError() const
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  std::ofstream csv(config_name_ + "-temperature.csv");
  csv << "x,y,temperature,exact,conduction_only,plane_conduction" << std::endl;
  double maximum{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double x{(*node_it)->x()};
    const double temperature{(*node_it)->Read(temperature_key)};
    maximum = std::max(maximum, std::fabs(temperature - ExactTemperature(x)));
    csv << x << "," << (*node_it)->y() << "," << temperature << "," << ExactTemperature(x) << ","
        << ConductionTemperature(x) << "," << PlaneTemperature(x) << std::endl;
  }
  return maximum;
}


void CVFEM_AxisymmetricAnnulus_Test::run()
{
  const double revolved_volume{CreateModel()};
  ReadGeometryAndWallTemperatures();
  const double exact_volume{M_PI * (outer_radius_ * outer_radius_ - inner_radius_ * inner_radius_) * height_};
  const double volume_error{std::fabs(revolved_volume - exact_volume) / exact_volume};

  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};
  InitializeAllNaNsToZero();
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(run_settings_.TimeIncrement() * year_);
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(run_settings_.TimeIncrement() * year_ / 10.);
  Advance(max_time);

  const std::array<double, 2> mass_rate{RadialMassRate()};
  const double heat_capacity_fluid{(WallEnthalpy(inner_radius_) - WallEnthalpy(outer_radius_))
                                   / (inner_temperature_ - outer_temperature_)};
  const double conductivity{(*model_->Region("Model").CellsBegin())->Read(
      property_database_->StorageKey("thermal conductivity"))};
  peclet_ = heat_capacity_fluid * mass_rate[0] / (2. * M_PI * height_ * conductivity);
  const double max_temperature_error{MaxTemperatureError()};

  // the element-centroid thickness integrates 1/r with the midpoint rule: about 0.6 K error next to the inner wall
  // of test-210a/b (1 K at the first structured node), against 34 K between the plane and the logarithmic profile
  const double temperature_limit{2.e-2 * std::fabs(inner_temperature_ - outer_temperature_)};
  double largest_conduction_difference{0.};
  for (double r = inner_radius_; r <= outer_radius_; r += 1.)
    largest_conduction_difference = std::max(largest_conduction_difference,
                                             std::fabs(ExactTemperature(r) - PlaneTemperature(r)));
  std::cerr << "CVFEM_AxisymmetricAnnulus_Test::run: " << mesh_name_ << ": r " << inner_radius_ << "-"
            << outer_radius_ << " m, height " << height_ << " m, " << steps_ << " steps to "
            << model_time / year_ << " yr" << std::endl
            << "  revolved volume " << revolved_volume << " m3, exact " << exact_volume << " m3 (relative error "
            << volume_error << ")" << std::endl
            << "  radial mass rate " << mass_rate[0] << " kg/s (largest element deviation of r rho q_r "
            << mass_rate[1] << "), fluid heat capacity " << heat_capacity_fluid << " J/kg/K, Peclet number A "
            << peclet_ << std::endl
            << "  largest |T - T_exact| " << max_temperature_error << " K (limit " << temperature_limit
            << " K); the plane conduction profile differs from T_exact by up to " << largest_conduction_difference
            << " K" << std::endl;

  do_test(model_time >= max_time, "Axisymmetric annulus: ran the configured duration", __FILE__, __LINE__);
  do_test(volume_error < 1.e-12, "Axisymmetric annulus: revolved bulk volume equals pi (r2^2 - r1^2) H",
          __FILE__, __LINE__);
  do_test(max_temperature_error < temperature_limit,
          "Axisymmetric annulus: steady temperature within 2 % of the wall difference of the exact profile",
          __FILE__, __LINE__);
}

} // namespace csmp
