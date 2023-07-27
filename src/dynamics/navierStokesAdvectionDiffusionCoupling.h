/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H

#include "core/operator.h"

namespace olb {


/// Coupling between a Navier-Stokes and an Advection-Diffusion lattice
struct NavierStokesAdvectionDiffusionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,T0>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }

};


/// AD coupling with Boussinesq bouancy for Smagorinsky-LES
struct SmagorinskyBoussinesqCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct PR_TURB : public descriptors::FIELD_BASE<1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,SMAGORINSKY_PREFACTOR,PR_TURB,T0,OMEGA_NSE,OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();

    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    auto u = cellADE.template getField<descriptors::VELOCITY>();
    // tau coupling
    auto tauNS = cellNSE.template getFieldPointer<descriptors::TAU_EFF>();
    auto tauAD = cellADE.template getFieldPointer<descriptors::TAU_EFF>();

    V rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cellNSE.computeAllMomenta(rho, u.data(), pi);
    cellADE.template setField<descriptors::VELOCITY>(u);
    V PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    V PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    V tau_mol_NS = V{1} / parameters.template get<OMEGA_NSE>();
    V tau_mol_AD = V{1} / parameters.template get<OMEGA_ADE>();
    /// Turbulent realaxation time
    V smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    V tau_turb_NS = V{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    /// Effective realaxation time
    tauNS[0] = tau_mol_NS+tau_turb_NS;

    V prTurb = parameters.template get<PR_TURB>();
    V tauTurbADPrefactor = descriptors::invCs2<V,DESCRIPTOR_ADE>() / descriptors::invCs2<V,DESCRIPTOR>() / prTurb;
    V tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
    tauAD[0] = tau_mol_AD+tau_turb_AD;
  }

};


}

#endif
