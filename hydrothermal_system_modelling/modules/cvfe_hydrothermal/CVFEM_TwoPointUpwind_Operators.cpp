// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_TwoPointUpwind_Operators.cpp
 *  @author Samuel Scott
 *  @brief Pressure-equation terms with two-point facet fluxes, see CVFEM_TwoPointUpwind_Operators.h.
 *  @date 20.11.2026
 */

#include "CVFEM_TwoPointUpwind_Operators.h"
#include "TwoPointTransmissibility.h"

namespace csmp {

namespace {

/// Mobility of the donor of the facet between inside and outside (0 without a donor).
double donorValue(const DenseMatrix<DM_MIN>& upwind, const std::vector<ScalarVariable>& values, uint32_t inside,
                  uint32_t outside)
{
  if (upwind(inside, outside) == 1) return values[inside]();
  if (upwind(inside, outside) == 2) return values[outside]();
  return 0.;
}

} // namespace

// ── conductance ──────────────────────────────────────────────────────────────

template<uint32_t dim>
CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>::CVFEM_TwoPointUpwind_dNT_op_dN_dV(
    const PropertyDatabase<dim>& pref, UpwindControlVisitor<dim>& upwind_visitor,
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport, const char* oper, const char* basic, const char* test,
    const char* upwind, const char* grav_trigger, const char* thickness)
  : CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element>(pref, upwind_visitor, fv_transport, oper, basic, test, upwind,
                                                        grav_trigger, thickness),
    upwind_visitor_(upwind_visitor),
    finite_volume_(fv_transport),
    mobility_(),
    mobility_multiplier_(),
    transmissibility_(),
    mobility_key_(pref.StorageKey(upwind)),
    gravity_trigger_key_(pref.StorageKey(grav_trigger)),
    thickness_key_(pref.StorageKey(thickness))
{
}


template<uint32_t dim>
CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>* CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>::clone() const
{
  return new CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>(*this);
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>::GetOperands(const Element<dim>& e)
{
  CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element>::GetOperands(e);
  e.NodePropertyVector(mobility_key_, mobility_);
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>::GetOperandsCVFEM(const Element<dim>& e, csmp::Index upwind_var_key)
{
  CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element>::GetOperandsCVFEM(e, upwind_var_key);
  e.NodePropertyVector(mobility_key_, mobility_);
  e.NodePropertyVector(upwind_var_key, mobility_multiplier_);
  for (size_t i{0U}; i < mobility_.size(); i++)
    mobility_[i]() *= mobility_multiplier_[i]();
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_dNT_op_dN_dV<dim>::ComputeContribution(const Element<dim>& e)
{
  TensorVariable<dim> permeability;   // diagonal, as the base operator
  permeability = 0.;
  for (uint32_t r{0U}; r < dim; r++)
    permeability(r, r) = this->MTRL[0](r, r);
  if (!finite_volume_.TwoPointFacetVelocities() || !twoPointTransmissibilities(e, permeability, transmissibility_)) {
    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim, Element>::ComputeContribution(e);
    return;
  }

  const DenseMatrix<DM_MIN>& upwind{upwind_visitor_.UpwindMatrix(gravity_trigger_key_, e.Idx())};
  DenseMatrix<DM_MIN>& lhs{this->LHS};
  lhs.Resize(e.Nodes(), e.Nodes());
  lhs.Zero();
  const double thickness{e.Read(thickness_key_)};
  for (uint32_t i{0U}; i < e.FV()->Facets(); i++) {
    uint32_t inside{0U}, outside{0U};
    e.FV()->FacetEdgeNodes(i, inside, outside);
    // flux out of inside = coefficient (p_inside - p_outside)
    const double coefficient{transmissibility_[i] * donorValue(upwind, mobility_, inside, outside) * thickness};
    lhs(inside, inside) += coefficient;
    lhs(inside, outside) -= coefficient;
    lhs(outside, inside) -= coefficient;
    lhs(outside, outside) += coefficient;
  }
}

// ── gravity ──────────────────────────────────────────────────────────────────

template<uint32_t dim>
CVFEM_TwoPointUpwind_rhsop_g_dV<dim>::CVFEM_TwoPointUpwind_rhsop_g_dV(
    const PropertyDatabase<dim>& pref, UpwindControlVisitor<dim>& upwind_visitor,
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport, const char* oper, const char* test, const char* upwind,
    const char* grav_trigger, const char* thickness)
  : CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>(pref, upwind_visitor, fv_transport, oper, test, upwind, grav_trigger,
                                                 thickness),
    upwind_visitor_(upwind_visitor),
    finite_volume_(fv_transport),
    mobility_density_(),
    mobility_multiplier_(),
    density_(),
    transmissibility_(),
    mobility_density_key_(pref.StorageKey(upwind)),
    density_key_(pref.StorageKey(grav_trigger)),
    thickness_key_(pref.StorageKey(thickness)),
    permeability_tensor_key_(pref.StorageKey("permeability tensor"))
{
}


template<uint32_t dim>
CVFEM_TwoPointUpwind_rhsop_g_dV<dim>* CVFEM_TwoPointUpwind_rhsop_g_dV<dim>::clone() const
{
  return new CVFEM_TwoPointUpwind_rhsop_g_dV<dim>(*this);
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_rhsop_g_dV<dim>::GetOperands(const Element<dim>& e)
{
  CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::GetOperands(e);
  e.NodePropertyVector(mobility_density_key_, mobility_density_);
  e.NodePropertyVector(density_key_, density_);
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_rhsop_g_dV<dim>::GetOperandsCVFEM(const Element<dim>& e, csmp::Index upwind_var_key)
{
  CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::GetOperandsCVFEM(e, upwind_var_key);
  e.NodePropertyVector(mobility_density_key_, mobility_density_);
  e.NodePropertyVector(density_key_, density_);
  e.NodePropertyVector(upwind_var_key, mobility_multiplier_);
  for (size_t i{0U}; i < mobility_density_.size(); i++)
    mobility_density_[i]() *= mobility_multiplier_[i]();
}


template<uint32_t dim>
void CVFEM_TwoPointUpwind_rhsop_g_dV<dim>::ComputeContribution(const Element<dim>& e)
{
  TensorVariable<dim> permeability;
  e.Read(permeability_tensor_key_, permeability);
  if (!finite_volume_.TwoPointFacetVelocities() || !twoPointTransmissibilities(e, permeability, transmissibility_)) {
    CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::ComputeContribution(e);
    return;
  }

  const DenseMatrix<DM_MIN>& upwind{upwind_visitor_.UpwindMatrix(density_key_, e.Idx())};
  std::vector<double>& rhs{this->RHS};
  rhs.assign(e.Nodes(), 0.);
  std::vector<double> g(dim, 0.);
  g[dim == 1U ? 0U : 1U] = -9.80665;   // as StencilProcessorPHX
  const double thickness{e.Read(thickness_key_)};
  for (uint32_t i{0U}; i < e.FV()->Facets(); i++) {
    uint32_t inside{0U}, outside{0U};
    e.FV()->FacetEdgeNodes(i, inside, outside);
    const double donor_density{donorValue(upwind, density_, inside, outside)};
    if (donor_density <= 0.)
      continue;
    // donor mobility x facet-mean density (the density of the two-point potential), as the two-point transport
    const double rho_i{density_[inside]()}, rho_o{density_[outside]()};
    const double rho{(rho_i > 0. && rho_o > 0.) ? 0.5 * (rho_i + rho_o) : donor_density};
    double height_term{0.};   // -g.(x_inside - x_outside)
    for (uint32_t r{0U}; r < dim; r++)
      height_term -= g[r] * (e.N(inside)->Coordinate()[r] - e.N(outside)->Coordinate()[r]);
    const double flux{transmissibility_[i] * thickness * donorValue(upwind, mobility_density_, inside, outside)
                      / donor_density * rho * height_term};   // inside -> outside
    rhs[inside] -= flux;
    rhs[outside] += flux;
  }
}


template class CVFEM_TwoPointUpwind_dNT_op_dN_dV<1U>;
template class CVFEM_TwoPointUpwind_dNT_op_dN_dV<2U>;
template class CVFEM_TwoPointUpwind_dNT_op_dN_dV<3U>;
template class CVFEM_TwoPointUpwind_rhsop_g_dV<1U>;
template class CVFEM_TwoPointUpwind_rhsop_g_dV<2U>;
template class CVFEM_TwoPointUpwind_rhsop_g_dV<3U>;

} // namespace csmp
