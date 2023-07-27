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

//This file contains the Local Pressure Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_LOCAL_PRESSURE_BOUNDARY_3D_HH
#define SET_LOCAL_PRESSURE_BOUNDARY_3D_HH

#include "setLocalPressureBoundary3D.h"

namespace olb {

///Initialising the setLocalPressureBoundary function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setLocalPressureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,3>& superGeometry, int material)
{
  setLocalPressureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega,superGeometry.getMaterialIndicator(material));
}

///Initialising the setLocalPressureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setLocalPressureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setLocalPressureBoundary");
  int _overlap = 0;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  //clout << sLattice.getLoadBalancer().size() <<"sLattice.getLoadBalancer.size()" << std::endl;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setLocalPressureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iC), omega, indicator->getBlockIndicatorF(iC),includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Set pressure boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setLocalPressureBoundary(BlockLattice<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

      if (discreteNormal[0] == 0) {//set momenta, dynamics and postProcessors on indicated cells
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 0,-1>>>();
        }
        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 0,1>>>();
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 1,-1>>>();
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 1,1>>>();
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 2,-1>>>();
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
            MixinDynamics,momenta::RegularizedPressureBoundaryTuple< 2,1>>>();
        }
      }
      if (dynamics) {
        dynamics->getParameters(_block).template set<descriptors::OMEGA>(omega);
      }
      setBoundary(_block, iX,iY,iZ, dynamics);
    }
  });
}


}
#endif
