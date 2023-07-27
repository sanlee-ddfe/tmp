/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

/* poiseuille2d.cpp:
 * This example examines a 2D Poseuille flow
 * It illustrates the computation of error norms.
 */

// the main code of the simulation is in poiseuille2d.h as it is also used by the
// example poiseuille2dEoc

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "poiseuille2d.h"

//Initialize Gnuplot
static Gnuplot<T> gplot("centerVelocity");

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );

  // Simulation parameter
  int N = 50;
  int status = readParameters(argc, argv, N, flowType, boundaryType);
  if (status != 0) {
    return status;
  }
  bool eoc = false;

  simulatePoiseuille(N, gplot, eoc);
}
