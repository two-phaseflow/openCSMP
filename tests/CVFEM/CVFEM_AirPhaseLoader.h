// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_AIR_PHASE_LOADER_H
#define CVFEM_AIR_PHASE_LOADER_H

/** @file CVFEM_AirPhaseLoader.h
 *  @author Samuel Scott
 *  @brief Puts air into a water-saturated CVFEM PHX model before the first time step (test-208, test-209).
 *  @details CVFEM_PHX_Scheme::InitialFluidPropertiesFromPTX sets every node to pure H2O-NaCl fluid with no air; in a
 *  model with an open top, air enters later through the boundary. The air-phase tests need air in the model from the
 *  start. Load() replaces part of the pore water by a given air mass and re-equilibrates every node once with a
 *  NaClH2OPropertiesVisitorPHX configured as the scheme's own (air in the hydraulic-mixture mode, thermodynamic
 *  density in the gravity and compressibility terms, fluid contribution to the heat capacity). The change is handed
 *  to the visitor as a transport increment: the liquid mass and enthalpy content fall to liquid_saturation times their
 *  water-saturated values and the air mass and enthalpy content rise from zero, so the visitor's own mass and
 *  enthalpy bookkeeping (total mass, total enthalpy, air saturation, mobilities, compressibility, pressure source)
 *  starts consistent. Pore water is pure water (no salt). The liquid saturation can vary by node (test-211).
 *  @date 20.11.2026
 */

#include "Model.h"
#include "NaClH2OPropertiesVisitorPHX.h"

#include <functional>

namespace csmp {

/// Replaces part of the pore water of a water-saturated CVFEM PHX model by air; total air mass of the model.
class CVFEM_AirPhaseLoader {
  public:
    /// model: initialised with CVFEM_PHX_Scheme::InitialFluidPropertiesFromPTX (air activated, water-saturated);
    /// rock_heat_capacity as passed to CVFEM_PHX_Scheme::SetRockHeatCapacity.
    CVFEM_AirPhaseLoader(Model<2U>& model, double rock_heat_capacity);
    ~CVFEM_AirPhaseLoader();

    /// Leaves liquid_saturation of the pore volume to the water, adds air_mass(node) [kg per m3 of pore volume] and
    /// re-equilibrates every node.
    void Load(double liquid_saturation, const std::function<double(const Node<2U>&)>& air_mass);
    /// As Load(double, ...) with a liquid saturation per node (test-211: water-saturated below the water table).
    void Load(const std::function<double(const Node<2U>&)>& liquid_saturation,
              const std::function<double(const Node<2U>&)>& air_mass);
    /// Air mass of the model [kg]: sum over the nodes of "fluid mass air" times "pore volume".
    double TotalAirMass() const;

  private:
    /// Liquid and air masses and enthalpy contents as a transport increment from the water-saturated state.
    void SetIncrementFromWaterSaturatedState(const std::function<double(const Node<2U>&)>& liquid_saturation,
                                             const std::function<double(const Node<2U>&)>& air_mass);

    Model<2U>& model_;
    NaClH2OPropertiesVisitorPHX<2U> properties_visitor_;
};

} // namespace csmp

#endif // CVFEM_AIR_PHASE_LOADER_H
