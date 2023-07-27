/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef SHAN_CHEN_FORCED_POST_PROCESSOR_H
#define SHAN_CHEN_FORCED_POST_PROCESSOR_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"

namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

struct RhoStatistics  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    auto statistic = cell.template getField<descriptors::STATISTIC>();
    statistic[0] = cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(statistic);
  }
};

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

template <typename POTENTIAL>
struct ShanChenForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct RHO0 : public descriptors::FIELD_BASE<2> { };
  struct G    : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_A : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_B : public descriptors::FIELD_BASE<1> { };

  using parameters = typename POTENTIAL::parameters::template decompose_into<
    meta::list<RHO0,G,OMEGA_A,OMEGA_B>::template include
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    auto& cellA = cells.template get<names::A>();
    auto& cellB = cells.template get<names::B>();

    auto rho0 = parameters.template get<RHO0>();

    Vector<V,2> rhoField{
      cellA.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[0],
      cellB.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[1]
    };
    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    V omegaA = parameters.template get<OMEGA_A>();
    V omegaB = parameters.template get<OMEGA_B>();
    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoField[0]*omegaA
             + rhoField[1]*omegaB;

    Vector<V,DESCRIPTOR::d> uTot = (uA*rho0[0]*omegaA + uB*rho0[1]*omegaB) / rhoTot;

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoBlockContribution{};
    Vector<V,DESCRIPTOR::d> rhoPartnerContribution{};
    V psi1 = POTENTIAL().compute(rhoField[0], parameters);
    V psi2 = POTENTIAL().compute(rhoField[1], parameters);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V rhoA = POTENTIAL().compute(nextCellA.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[0], parameters);
      V rhoB = POTENTIAL().compute(nextCellB.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[1], parameters);
      rhoBlockContribution   += psi2 * rhoA * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      rhoPartnerContribution += psi1 * rhoB * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    auto externalBlockForce   = cellA.template getField<descriptors::EXTERNAL_FORCE>();
    auto externalPartnerForce = cellB.template getField<descriptors::EXTERNAL_FORCE>();

    auto g = parameters.template get<G>();

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);
    cellA.template setField<descriptors::FORCE>(externalBlockForce
        - g*rhoPartnerContribution/rhoField[0]);
    cellB.template setField<descriptors::FORCE>(externalPartnerForce
        - g*rhoBlockContribution/rhoField[1]);
  }

};

}

#endif
