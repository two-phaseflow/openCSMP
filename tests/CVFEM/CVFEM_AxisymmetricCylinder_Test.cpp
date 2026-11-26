// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_AxisymmetricCylinder_Test.cpp
 *  @author Samuel Scott
 *  @brief A solid cylinder cooling by conduction, axis inside the mesh, see CVFEM_AxisymmetricCylinder_Test.h.
 *  @date 26.11.2026
 */

#include "CVFEM_AxisymmetricCylinder_Test.h"

#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "PoreVolumeVisitor.h"
#include "WellConfigurationFile.h"
#include "ModelTime.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace csmp {

CVFEM_AxisymmetricCylinder_Test::CVFEM_AxisymmetricCylinder_Test(const std::string& mesh_name,
                                                                 const std::string& config_name)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    steps_(0U),
    radius_(0.),
    height_(0.),
    start_temperature_(0.),
    wall_temperature_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    run_settings_(),
    axisymmetric_{true, 0., M_PI}   // axis inside the mesh: each half sweeps 180 degrees
{
}


CVFEM_AxisymmetricCylinder_Test::~CVFEM_AxisymmetricCylinder_Test() = default;


double CVFEM_AxisymmetricCylinder_Test::CreateModel()
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


void CVFEM_AxisymmetricCylinder_Test::InitializeAllNaNsToZero()
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


void CVFEM_AxisymmetricCylinder_Test::ReadGeometryAndTemperatures()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  double y_bottom{std::numeric_limits<double>::max()};
  double y_top{std::numeric_limits<double>::lowest()};
  start_temperature_ = std::numeric_limits<double>::lowest();
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    radius_ = std::max(radius_, std::fabs((*node_it)->x() - axisymmetric_.axis_x));
    y_bottom = std::min(y_bottom, (*node_it)->y());
    y_top = std::max(y_top, (*node_it)->y());
    start_temperature_ = std::max(start_temperature_, (*node_it)->Read(temperature_key));
  }
  height_ = y_top - y_bottom;
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if (std::fabs(std::fabs((*node_it)->x() - axisymmetric_.axis_x) - radius_) < 1.e-6)
      wall_temperature_ = (*node_it)->Read(temperature_key);
}


void CVFEM_AxisymmetricCylinder_Test::Advance(double max_time)
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


double CVFEM_AxisymmetricCylinder_Test::Diffusivity() const
{
  const Index porosity_key(property_database_->StorageKey("nodal porosity"));
  const Index rock_density_key(property_database_->StorageKey("nodal density rock"));
  const Index rock_heat_capacity_key(property_database_->StorageKey("nodal heat capacity rock"));
  const Index fluid_density_key(property_database_->StorageKey("density liquid"));
  const Index fluid_heat_capacity_key(property_database_->StorageKey("fluid heat capacity"));
  const Region<2U>& region(model_->Region("Model"));
  double heat_capacity{0.};
  uint32_t count{0U};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double porosity{(*node_it)->Read(porosity_key)};
    heat_capacity += (1. - porosity) * (*node_it)->Read(rock_density_key) * (*node_it)->Read(rock_heat_capacity_key)
                     + porosity * (*node_it)->Read(fluid_density_key) * (*node_it)->Read(fluid_heat_capacity_key);
    ++count;
  }
  const double conductivity{(*model_->Region("Model").CellsBegin())->Read(
      property_database_->StorageKey("thermal conductivity"))};
  return conductivity / (heat_capacity / count);
}


double CVFEM_AxisymmetricCylinder_Test::ExactTemperature(double r, double fourier) const
{
  // first zeros of J0; the omitted terms are below exp(-17.8^2 * 0.05) ~ 1e-7 for Fo >= 0.05
  constexpr std::array<double, 6> zeros{2.404825557695773, 5.520078110286311, 8.653727912911013,
                                        11.79153443901428, 14.93091770848779, 18.07106396791092};
  double sum{0.};
  for (double zero : zeros)
    sum += 2. / (zero * std::cyl_bessel_j(1., zero)) * std::cyl_bessel_j(0., zero * r / radius_)
           * std::exp(-zero * zero * fourier);
  return wall_temperature_ + (start_temperature_ - wall_temperature_) * sum;
}


double CVFEM_AxisymmetricCylinder_Test::PlaneTemperature(double x, double fourier) const
{
  double sum{0.};
  for (int n = 0; n < 12; ++n) {
    const double k{(2. * n + 1.) * M_PI / 2.};
    sum += 2. * (n % 2 == 0 ? 1. : -1.) / k * std::cos(k * x / radius_) * std::exp(-k * k * fourier);
  }
  return wall_temperature_ + (start_temperature_ - wall_temperature_) * sum;
}


void CVFEM_AxisymmetricCylinder_Test::run()
{
  const double revolved_volume{CreateModel()};

  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};
  InitializeAllNaNsToZero();
  ReadGeometryAndTemperatures();
  const double exact_volume{M_PI * radius_ * radius_ * height_};
  const double volume_error{std::fabs(revolved_volume - exact_volume) / exact_volume};

  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(run_settings_.TimeIncrement() * year_);
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(run_settings_.TimeIncrement() * year_ / 10.);
  Advance(max_time);

  const double fourier{Diffusivity() * model_time / (radius_ * radius_)};
  const Index temperature_key(property_database_->StorageKey("temperature"));
  std::ofstream csv(config_name_ + "-temperature.csv");
  csv << "x,y,temperature,exact,plane" << std::endl;
  double max_error{0.};
  double max_plane_difference{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double x{(*node_it)->x() - axisymmetric_.axis_x};
    const double temperature{(*node_it)->Read(temperature_key)};
    const double exact{ExactTemperature(std::fabs(x), fourier)};
    max_error = std::max(max_error, std::fabs(temperature - exact));
    max_plane_difference = std::max(max_plane_difference, std::fabs(PlaneTemperature(x, fourier) - exact));
    csv << (*node_it)->x() << "," << (*node_it)->y() << "," << temperature << "," << exact << ","
        << PlaneTemperature(x, fourier) << std::endl;
  }
  const double temperature_limit{2.e-2 * std::fabs(start_temperature_ - wall_temperature_)};
  std::cerr << "CVFEM_AxisymmetricCylinder_Test::run: " << mesh_name_ << ": R " << radius_ << " m, height " << height_
            << " m, " << steps_ << " steps to " << model_time / year_ << " yr (Fourier number " << fourier << ")"
            << std::endl
            << "  revolved volume " << revolved_volume << " m3, exact " << exact_volume << " m3 (relative error "
            << volume_error << ")" << std::endl
            << "  largest |T - T_exact| " << max_error << " K (limit " << temperature_limit
            << " K); the plane solution differs from T_exact by up to " << max_plane_difference << " K" << std::endl;

  do_test(model_time >= max_time, "Axisymmetric cylinder: ran the configured duration", __FILE__, __LINE__);
  do_test(volume_error < 1.e-12, "Axisymmetric cylinder: revolved bulk volume equals pi R^2 H", __FILE__, __LINE__);
  do_test(max_error < temperature_limit,
          "Axisymmetric cylinder: temperature within 2 % of T0 - T2 of the Bessel series", __FILE__, __LINE__);
}

} // namespace csmp
