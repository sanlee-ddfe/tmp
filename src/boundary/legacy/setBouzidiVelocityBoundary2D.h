/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

//This file contains the Bouzidi Velocity Boundary
//This is an offLattice boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_BOUZIDI_VELOCITY_BOUNDARY_2D_H
#define SET_BOUZIDI_VELOCITY_BOUNDARY_2D_H

#include <vector>
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry.h"
#include "utilities/functorPtr.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/freeEnergyDynamics.h"
#include "boundary/offBoundaryPostProcessors2D.h"

namespace olb {

namespace legacy {

////////// SuperLattice Domain  /////////////////////////////////////////

/// Set offDynamics with boundary links and post processors using indicators
/**
 * Add offDynamics with initialisation of boundary links and the corresponding
 * post processors
 * Note: Uses information of the second neighbours of the cell (x,y)
 * Add post processors. Ensure that offDynamics are defined!
 *
 * \param boundaryIndicator Indicator describing boundary cells
 * \param bulkIndicator     Indicator describing bulk cells
 * \param geometryIndicator Indicator describing the geometry to be bounded
 **/

///Initialising the BouzidiVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics = BGKdynamics<T,DESCRIPTOR>>
void setBouzidiVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material,
                                IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));

///Initialising the BouzidiVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                                FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                                IndicatorF2D<T>& geometryIndicator);

////////// BlockLattice Domain  /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator);

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockGeometry<T,2>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator);

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockGeometry<T,2>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q]);

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockGeometry<T,2>& blockGeometryStructure, int x, int y, int iPop, T dist);

//set dynamics on indicated cells
//this function can be used for bouzidi and bounceBackVelocityBoundary
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setOffDynamics(BlockLattice<T,DESCRIPTOR>& block, int x, int y, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]);

}

}//namespace olb

#endif
