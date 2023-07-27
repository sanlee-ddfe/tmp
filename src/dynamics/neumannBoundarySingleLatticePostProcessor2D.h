/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Johanna Moedl
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

#ifndef NEUMANN_BOUNDARY_SINGLE_LATTICE_POST_PROCESSOR_2D_H
#define NEUMANN_BOUNDARY_SINGLE_LATTICE_POST_PROCESSOR_2D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "latticeDescriptors.h"
#include "utilities/omath.h"


namespace olb {


//======================================================================
// ======== Neumann Boundary with Finite Differences for AD 2D ======//
//======================================================================
template<typename T, typename DESCRIPTOR, int normal1, int normal2>
struct NeumannBoundarySingleLatticePostProcessor2D {
 static constexpr OperatorScope scope = OperatorScope::PerCell;

 int getPriority() const {
   return 0;
 }
 template <typename CELL>
 void apply(CELL& cell) any_platform {
   T neumannBoundary;
   T analyticalBoundaryValue = cell.template getField<descriptors::BOUNDARY>();

   if constexpr (normal2 == 0 && normal1 != 0) {
     if constexpr (normal1 > 0) {
       //right boundary, difference quotient with x-1
       neumannBoundary = analyticalBoundaryValue+cell.neighbor({-1,0}).computeRho();
     }
     else if constexpr (normal1<0){
       //left boundary, difference quotient with x+1
       neumannBoundary = cell.neighbor({1,0}).computeRho() - analyticalBoundaryValue;
     }
   }

   if constexpr (normal1 == 0 && normal2 != 0) {
     if constexpr (normal2 > 0) {
       //upper boundary, difference quotient with y-1
       neumannBoundary = analyticalBoundaryValue + cell.neighbor({0,-1}).computeRho();
     }
     else if constexpr (normal2<0){
       //lower boundary, difference quotient with y+1
       neumannBoundary = cell.neighbor({0,1}).computeRho() + analyticalBoundaryValue;
     }
   }

   cell.defineRho(neumannBoundary);
 }
};

}

#endif
