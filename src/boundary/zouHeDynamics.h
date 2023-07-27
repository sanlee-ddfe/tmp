/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef ZOU_HE_DYNAMICS_H
#define ZOU_HE_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of Zou-He boundary condition following
* the paper from Zou and He. This implementation is lattice independent.
* The implementation follow the idea proposed int the paper
* Qisu Zou, Xiaoyi He,
* "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model",
* Phys. Fluids , (1997), Volume 9, Issue 6, pp. 1591-1598
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
class ZouHeDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template<typename M>
  using exchange_momenta = ZouHeDynamics<T,DESCRIPTOR,DYNAMICS,M,direction,orientation>;

  std::type_index id() override {
    return typeid(ZouHeDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ZouHeDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    // Along all the commented parts of this code there will be an example based
    // on the situation where the wall's normal vector if (0,1) and the
    // numerotation of the velocites are done according to the D2Q9
    // lattice of the OpenLB library.

    // Find all the missing populations
    // (directions 3,4,5)
    constexpr auto missingIndexesTmp = util::subIndexOutgoing<DESCRIPTOR,direction,orientation>();
    std::vector<int> missingIndexes(missingIndexesTmp.cbegin(), missingIndexesTmp.cend());

    // Will contain the missing poputations that are not normal to the wall.
    // (directions 3,5)
    std::vector<int> missingDiagonalIndexes = missingIndexes;
    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
      int numOfNonNullComp = 0;
      for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
        numOfNonNullComp += util::abs(descriptors::c<DESCRIPTOR>(missingIndexes[iPop],iDim));
      }

      if (numOfNonNullComp == 1) {
        missingDiagonalIndexes.erase(missingDiagonalIndexes.begin()+iPop);
        break;
      }
    }

    V rho, u[DESCRIPTOR::d];
    V falseRho, falseU[DESCRIPTOR::d];
    MomentaF().computeRhoU(cell, rho, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

    // The unknown non equilibrium populations are bounced back
    // (f[3] = feq[3] + fneq[7], f[4] = feq[4] + fneq[8],
    //  f[5] = feq[5] + fneq[1])
    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
      cell[missingIndexes[iPop]] = cell[descriptors::opposite<DESCRIPTOR>(missingIndexes[iPop])]
        - computeEquilibrium(descriptors::opposite<DESCRIPTOR>(missingIndexes[iPop]), rho, u)
        + computeEquilibrium(missingIndexes[iPop], rho, u);
    }

    // We recompute rho and u in order to have the new momentum and density. Since
    // the momentum is not conserved from this scheme, we will corect it. By adding
    // a contribution to the missingDiagonalVelocities.
    lbm<DESCRIPTOR>::computeRhoU(cell,falseRho,falseU);

    V diff[DESCRIPTOR::d] {};
    for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
      diff[iDim] = (rho*u[iDim] - falseRho*falseU[iDim])/ V(missingDiagonalIndexes.size());
    }

    for (unsigned iPop = 0; iPop < missingDiagonalIndexes.size(); ++iPop) {
      for (int iDim = 1; iDim < DESCRIPTOR::d; ++iDim) {
        cell[missingDiagonalIndexes[iPop]] +=
          descriptors::c<DESCRIPTOR>(missingDiagonalIndexes[iPop],(direction+iDim)%DESCRIPTOR::d) * diff[(direction+iDim)%DESCRIPTOR::d];
      }
    }

    typename CORRECTED_DYNAMICS::CollisionO().apply(cell, parameters);

    return {rho, uSqr};
  }

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override any_platform {
    return typename CORRECTED_DYNAMICS::EquilibriumF().compute(iPop, rho, u);
  };

  std::string getName() const override {
    return "ZouHeDynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

}

#endif
