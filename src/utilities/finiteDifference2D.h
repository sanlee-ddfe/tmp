/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef FINITE_DIFFERENCE_2D_H
#define FINITE_DIFFERENCE_2D_H

#include "finiteDifference.h"

namespace olb {

namespace fd {

// -------------------- Second derivatives ------------------------------------
// cf. Abramowitz, Stegun: Handbook of Mathematical Functions, 1964, p. 884f.

/// Second order mixed derivative (u_xp1_ym1 = u(x+1,y-1))
template<typename T>
constexpr T d2u_dxdy(T u_xm1_ym1, T u_xm1_yp1, T u_xp1_ym1, T u_xp1_yp1) any_platform
{
  return (u_xp1_yp1 - u_xp1_ym1 - u_xm1_yp1 + u_xm1_ym1) / T{4};
}

/// Second order mixed derivative (u_xp1_ym1 = u(x+1,y-1))
template<typename T>
constexpr T d2u_dxdy(T u_xm1_ym1, T u_xm1, T u_ym1,
  T u_0, T u_xp1, T u_yp1, T u_xp1_yp1) any_platform
{
  return - (u_xp1 + u_yp1 + u_xm1 + u_ym1 - T{2}*u_0 - u_xp1_yp1 - u_xm1_ym1) / T{2};
}

/// Second order Laplacian (u_xp1 = u(x+1,y))
template<typename T>
constexpr T laplacian2D(T u_xm1, T u_ym1, T u_0, T u_xp1, T u_yp1) any_platform
{
  return u_xm1 + u_ym1 + u_xp1 + u_yp1 - T{4}*u_0;
}

/// Forth order Laplacian (u_xp1 = u(x+1,y))
template<typename T>
constexpr T laplacian2D(T u_xm2, T u_ym2, T u_xm1, T u_ym1,
  T u_0, T u_xp1, T u_yp1, T u_xp2, T u_yp2) any_platform
{
  return (-T{60}*u_0 + T{16}*(u_xm1 + u_ym1 + u_xp1 + u_yp1)
         - (u_xm2 + u_ym2 + u_xp2 + u_yp2)) / T{12};
}


// -------------------- Interpolation -----------------------------------------

template<typename T, typename DESCRIPTOR,
         int direction, int orientation,
         bool orthogonal>
struct DirectedGradients2D {
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY);
  static void interpolateScalar(T& rhoDeriv,
                                BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY);

  template <typename CELL>
  static void interpolateVector( T velDeriv[DESCRIPTOR::d],
                                 CELL& cell ) any_platform;
};

// Implementation for orthogonal==true; i.e. the derivative is along
// the boundary normal.
template<typename T, typename DESCRIPTOR,
         int direction, int orientation>
struct DirectedGradients2D<T, DESCRIPTOR, direction, orientation, true> {
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY)
  {
    using namespace fd;

    T u0[DESCRIPTOR::d], u1[DESCRIPTOR::d], u2[DESCRIPTOR::d];

    blockLattice.get(iX,iY).computeU(u0);
    blockLattice.get (
      iX+(direction==0 ? (-orientation):0),
      iY+(direction==1 ? (-orientation):0) ).computeU(u1);
    blockLattice.get (
      iX+(direction==0 ? (-2*orientation):0),
      iY+(direction==1 ? (-2*orientation):0) ).computeU(u2);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = -orientation * boundaryGradient(u0[iD], u1[iD], u2[iD]);
    }
  }

  static void interpolateScalar(T& rhoDeriv,
                                BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY)
  {
    using namespace fd;

    T rho0 = blockLattice.get(iX,iY).computeRho();
    T rho1 = blockLattice.get (
               iX+(direction==0 ? (-orientation):0),
               iY+(direction==1 ? (-orientation):0) ).computeRho();
    T rho2 = blockLattice.get (
               iX+(direction==0 ? (-2*orientation):0),
               iY+(direction==1 ? (-2*orientation):0) ).computeRho();

    rhoDeriv = -orientation * boundaryGradient(rho0, rho1, rho2);

  }

  template <typename CELL>
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                CELL& cell) any_platform
  {
    using namespace fd;

    T u0[DESCRIPTOR::d], u1[DESCRIPTOR::d], u2[DESCRIPTOR::d];

    cell.computeU(u0);
    cell.neighbor({(direction==0 ? (-orientation):0),
                   (direction==1 ? (-orientation):0)}).computeU(u1);
    cell.neighbor({(direction==0 ? (-2*orientation):0),
                   (direction==1 ? (-2*orientation):0)}).computeU(u2);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = -orientation * boundaryGradient(u0[iD], u1[iD], u2[iD]);
    }
  }
};


// Implementation for orthogonal==false; i.e. the derivative is aligned
// with the boundary.
template<typename T, typename DESCRIPTOR,
         int direction, int orientation>
struct DirectedGradients2D<T, DESCRIPTOR, direction, orientation, false> {
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY)
  {
    using namespace fd;

    T u_p1[DESCRIPTOR::d], u_m1[DESCRIPTOR::d];

    int deriveDirection = 1-direction;
    blockLattice.get (
      iX+(deriveDirection==0 ? 1:0),
      iY+(deriveDirection==1 ? 1:0) ).computeU(u_p1);
    blockLattice.get (
      iX+(deriveDirection==0 ? (-1):0),
      iY+(deriveDirection==1 ? (-1):0) ).computeU(u_m1);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = fd::centralGradient(u_p1[iD],u_m1[iD]);
    }
  }

  static void  interpolateScalar(T& rhoDeriv,
                                 BlockLattice<T,DESCRIPTOR> const& blockLattice,
                                 int iX, int iY)
  {
    using namespace fd;

    int deriveDirection = 1-direction;
    T rho_p1 = blockLattice.get (
                 iX+(deriveDirection==0 ? 1:0),
                 iY+(deriveDirection==1 ? 1:0) ).computeRho();
    T rho_m1 = blockLattice.get (
                 iX+(deriveDirection==0 ? (-1):0),
                 iY+(deriveDirection==1 ? (-1):0) ).computeRho();

    rhoDeriv = centralGradient(rho_p1, rho_m1);

  }

  template <typename CELL>
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                CELL& cell) any_platform
  {
    using namespace fd;

    T u_p1[DESCRIPTOR::d], u_m1[DESCRIPTOR::d];

    int deriveDirection = 1-direction;
    cell.neighbor({(deriveDirection==0 ? 1:0),
                   (deriveDirection==1 ? 1:0)}).computeU(u_p1);
    cell.neighbor({(deriveDirection==0 ? (-1):0),
                   (deriveDirection==1 ? (-1):0)}).computeU(u_m1);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = centralGradient(u_p1[iD],u_m1[iD]);
    }
  }
};

}  // namespace fd

}  // namespace olb


#endif
