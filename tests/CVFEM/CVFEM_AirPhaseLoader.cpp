// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_AirPhaseLoader.cpp
 *  @author Samuel Scott
 *  @brief Puts air into a water-saturated CVFEM PHX model, see CVFEM_AirPhaseLoader.h.
 *  @date 20.11.2026
 */

#include "CVFEM_AirPhaseLoader.h"

#include "Air.h"

namespace csmp {

CVFEM_AirPhaseLoader::CVFEM_AirPhaseLoader(Model<2U>& model, double rock_heat_capacity)
  : model_(model),
    properties_visitor_(model)
{
  // as CVFEM_PHX_Scheme with thermodynamic_density_in_gravity_term = true, SetRockHeatCapacity,
  // AddFluidContributionToHeatCapacity(true) and ActivateAir(0., false, true, 0.)
  properties_visitor_.SetThermodynamicDensityInGravityTermTo(true);
  properties_visitor_.SetThermodynamicDensityInCompressibilityTermTo(true);
  properties_visitor_.SetRockHeatCapacity(rock_heat_capacity);
  properties_visitor_.AddFluidContributionToHeatCapacity(true);
  properties_visitor_.ActivateAir(0., false, true, 0.);
}


CVFEM_AirPhaseLoader::~CVFEM_AirPhaseLoader() = default;


void CVFEM_AirPhaseLoader::SetIncrementFromWaterSaturatedState(
    const std::function<double(const Node<2U>&)>& liquid_saturation,
    const std::function<double(const Node<2U>&)>& air_mass)
{
  const PropertyDatabase<2U>& database(model_.Database());
  const Index fluid_mass_key(database.StorageKey("fluid density"));
  const Index liquid_specific_enthalpy_key(database.StorageKey("enthalpy liquid"));
  const Index liquid_mass_key(database.StorageKey("fluid mass liquid"));
  const Index previous_liquid_mass_key(database.StorageKey("previous fluid mass liquid"));
  const Index liquid_enthalpy_key(database.StorageKey("enthalpy content liquid"));
  const Index previous_liquid_enthalpy_key(database.StorageKey("previous enthalpy content liquid"));
  const Index air_mass_key(database.StorageKey("fluid mass air"));
  const Index previous_air_mass_key(database.StorageKey("previous fluid mass air"));
  const Index air_enthalpy_key(database.StorageKey("enthalpy content air"));
  const Index previous_air_enthalpy_key(database.StorageKey("previous enthalpy content air"));
  const Index temperature_key(database.StorageKey("temperature"));
  Air air;
  const Region<2U>& region(model_.Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    Node<2U>& node(**node_it);
    const auto store_value = [&node](const Index& key, double value) {   // keeps the flag
      ScalarVariable variable;
      node.Read(key, variable);
      variable() = value;
      node.Store(key, variable);
    };
    // water-saturated: all fluid mass is liquid ("fluid mass liquid" is not set by InitialFluidPropertiesFromPTX)
    const double liquid_mass{node.Read(fluid_mass_key)};
    const double liquid_enthalpy{liquid_mass * node.Read(liquid_specific_enthalpy_key)};
    const double mass_of_air{air_mass(node)};
    const double saturation{liquid_saturation(node)};
    store_value(previous_liquid_mass_key, liquid_mass);
    store_value(liquid_mass_key, saturation * liquid_mass);
    store_value(previous_liquid_enthalpy_key, liquid_enthalpy);
    store_value(liquid_enthalpy_key, saturation * liquid_enthalpy);
    store_value(previous_air_mass_key, 0.);
    store_value(air_mass_key, mass_of_air);
    store_value(previous_air_enthalpy_key, 0.);
    store_value(air_enthalpy_key, mass_of_air * air.Enthalpy(node.Read(temperature_key)));
  }
}


void CVFEM_AirPhaseLoader::Load(double liquid_saturation, const std::function<double(const Node<2U>&)>& air_mass)
{
  Load([liquid_saturation](const Node<2U>&) { return liquid_saturation; }, air_mass);
}


void CVFEM_AirPhaseLoader::Load(const std::function<double(const Node<2U>&)>& liquid_saturation,
                                const std::function<double(const Node<2U>&)>& air_mass)
{
  // no conduction increment: the visitor credits rock (and fluid) heat for temperature - previous temperature
  model_.CopyReplace("temperature", "previous temperature");
  SetIncrementFromWaterSaturatedState(liquid_saturation, air_mass);
  model_.Accept(properties_visitor_);
}


double CVFEM_AirPhaseLoader::TotalAirMass() const
{
  const Index air_mass_key(model_.Database().StorageKey("fluid mass air"));
  const Index pore_volume_key(model_.Database().StorageKey("pore volume"));
  double total{0.};
  const Region<2U>& region(model_.Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    total += (*node_it)->Read(air_mass_key) * (*node_it)->Read(pore_volume_key);
  return total;
}

} // namespace csmp
