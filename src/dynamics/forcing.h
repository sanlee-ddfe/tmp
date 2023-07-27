/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_FORCING_H
#define DYNAMICS_FORCING_H

#include "lbm.h"
#include "descriptorField.h"
#include "core/latticeStatistics.h"

namespace olb {

/// Dynamics combination rules for various forcing schemes
namespace forcing {

/// Dynamics combination rule implementing the forcing scheme by Guo et al.
template <template <typename> typename Forced = momenta::Forced>
struct Guo {
  static std::string getName() {
    return "GuoForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Guo et al.
//template <template <typename> typename Sourced = momenta::Sourced>
struct AdeGuo {
  static std::string getName() {
    return "AdeGuoForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      const auto u = cell.template getField<descriptors::VELOCITY>();
      const V rho = MomentaF().computeRho(cell);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto source = cell.template getField<descriptors::SOURCE>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V sourceTerm{};
      sourceTerm = source;
      sourceTerm *= descriptors::t<V,DESCRIPTOR>(iPop);
      sourceTerm *= (V{1} - omega * V{0.5});
      cell[iPop] += sourceTerm;
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Kupershtokh et al.
struct Kupershtokh {
  static std::string getName() {
    return "KupershtokhForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V u[DESCRIPTOR::d];
      MomentaF().computeU(cell, u);
      const auto statistic = CollisionO().apply(cell, parameters);
      const auto force = cell.template getField<descriptors::FORCE>();
      V uPlusDeltaU[DESCRIPTOR::d];
      for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
        uPlusDeltaU[iVel] = u[iVel] + force[iVel];
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += EquilibriumF().compute(iPop, statistic.rho, uPlusDeltaU)
                    - EquilibriumF().compute(iPop, statistic.rho, u);
      }
      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Shan and Chen
class ShanChen {
private:
  template <typename EQUILIBRIUM>
  struct VelocityShiftedEquilibrium {
    static std::string getName() {
      return "ShanChen<" + EQUILIBRIUM::getName() + ">";
    }

    template <typename DESCRIPTOR, typename MOMENTA>
    struct type {
      using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
      using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

      template <typename RHO, typename U>
      auto compute(int iPop, const RHO& rho, const U& u) any_platform {
        return EquilibriumF().compute(iPop, rho, u);
      }

      template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
      CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
        V rho, u[DESCRIPTOR::d], uSqr{};
        MomentaF().computeRhoU(cell, rho, u);
        const auto force = cell.template getFieldPointer<descriptors::FORCE>();
        for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
          u[iVel] += force[iVel] / parameters.template get<descriptors::OMEGA>();
          uSqr += u[iVel] * u[iVel];
        }
        for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          fEq[iPop] = EquilibriumF().compute(iPop, rho, u);
        }
        return {rho, uSqr};
      };
    };
  };

public:
  static std::string getName() {
    return "ShanChenForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename VelocityShiftedEquilibrium<EQUILIBRIUM>::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = combined_momenta<DESCRIPTOR,MOMENTA>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,VelocityShiftedEquilibrium<EQUILIBRIUM>>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d] { };
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

struct LinearVelocity {
  static std::string getName() {
    return "LinearVelocityForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
      MomentaF().computeAllMomenta(cell, rho, u, pi);
      auto force = cell.template getFieldPointer<descriptors::FORCE>();
      int nDim = DESCRIPTOR::d;
      V forceSave[nDim];
      // adds a+Bu to force, where
      //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
      //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
      auto v = cell.template getFieldPointer<descriptors::V12>();
      for (int iDim=0; iDim<nDim; ++iDim) {
        forceSave[iDim] = force[iDim];
        force[iDim] += v[iDim];
        for (int jDim=0; jDim<nDim; ++jDim) {
          force[iDim] += v[jDim + iDim*nDim + nDim]*u[jDim];
        }
      }
      for (int iVel=0; iVel<nDim; ++iVel) {
        u[iVel] += force[iVel] / V{2.};
      }

      auto statistics = CollisionO().apply(cell, parameters);
      V newOmega = parameters.template get<descriptors::OMEGA>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, newOmega, force);
      // Writing back to froce fector
      for (int iVel=0; iVel<nDim; ++iVel) {
        force[iVel] = forceSave[iVel];
      }
      return statistics;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

}

#endif
