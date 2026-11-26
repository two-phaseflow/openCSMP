// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_TWO_POINT_UPWIND_OPERATORS_H
#define CVFEM_TWO_POINT_UPWIND_OPERATORS_H

/** @file CVFEM_TwoPointUpwind_Operators.h
 *  @author Samuel Scott
 *  @brief Pressure-equation conductance and gravity terms of a phase whose transport uses two-point facet fluxes.
 *  @details Drop-in replacements for CVFEM_Upwind_NumIntegral_dNT_op_dN_dV (conductance) and
 *  CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV (gravity): same constructor arguments, same behaviour, until the phase's
 *  ExplicitFiniteVolumeTransportPHX is switched to two-point facet velocities (UseTwoPointFacetVelocities). Then, on
 *  linear simplices, the flux through the facet of edge (a, b) is T_ab lambda (Phi_a - Phi_b), Phi = p - rho g.x, with
 *  T_ab from TwoPointTransmissibility.h, lambda the donor's mobility (the same upwind map as the base operators) and
 *  rho the mean density of the two nodes. The conductance puts T_ab lambda (p_a - p_b) into the matrix, the gravity
 *  term T_ab lambda rho (-g.(x_a - x_b)) into the right-hand side, so the pressure solve and the transport use the
 *  same fluxes. Other elements fall back to the base operators.
 *  @date 20.11.2026
 */

#include "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h"
#include "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h"
#include "ExplicitFiniteVolumeTransportPHX.h"

#include <vector>

namespace csmp {

/// Conductance term (matrix); arguments as CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.
template<uint32_t dim>
class CVFEM_TwoPointUpwind_dNT_op_dN_dV : public CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element> {
  public:
    CVFEM_TwoPointUpwind_dNT_op_dN_dV(const PropertyDatabase<dim>& pref, UpwindControlVisitor<dim>& upwind_visitor,
                                      ExplicitFiniteVolumeTransportPHX<dim>& fv_transport, const char* oper,
                                      const char* basic, const char* test, const char* upwind,
                                      const char* grav_trigger, const char* thickness);

    void GetOperands(const Element<dim>& e) override;
    void GetOperandsCVFEM(const Element<dim>& e, csmp::Index upwind_var_key) override;
    void ComputeContribution(const Element<dim>& e) override;
    CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>* clone() const override;

  private:
    UpwindControlVisitor<dim>& upwind_visitor_;
    ExplicitFiniteVolumeTransportPHX<dim>& finite_volume_;
    std::vector<ScalarVariable> mobility_;
    std::vector<ScalarVariable> mobility_multiplier_;
    std::vector<double> transmissibility_;
    csmp::Index mobility_key_, gravity_trigger_key_, thickness_key_;
};

/// Gravity term (right-hand side); arguments as CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.
template<uint32_t dim>
class CVFEM_TwoPointUpwind_rhsop_g_dV : public CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim> {
  public:
    CVFEM_TwoPointUpwind_rhsop_g_dV(const PropertyDatabase<dim>& pref, UpwindControlVisitor<dim>& upwind_visitor,
                                    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport, const char* oper,
                                    const char* test, const char* upwind, const char* grav_trigger,
                                    const char* thickness);

    void GetOperands(const Element<dim>& e) override;
    void GetOperandsCVFEM(const Element<dim>& e, csmp::Index upwind_var_key) override;
    void ComputeContribution(const Element<dim>& e) override;
    CVFEM_TwoPointUpwind_rhsop_g_dV<dim>* clone() const override;

  private:
    UpwindControlVisitor<dim>& upwind_visitor_;
    ExplicitFiniteVolumeTransportPHX<dim>& finite_volume_;
    std::vector<ScalarVariable> mobility_density_;   ///< donor "mass mobility density" (mobility x density)
    std::vector<ScalarVariable> mobility_multiplier_;
    std::vector<ScalarVariable> density_;
    std::vector<double> transmissibility_;
    csmp::Index mobility_density_key_, density_key_, thickness_key_, permeability_tensor_key_;
};

} // namespace csmp

#endif // CVFEM_TWO_POINT_UPWIND_OPERATORS_H
