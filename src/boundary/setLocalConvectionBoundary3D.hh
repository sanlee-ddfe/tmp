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

//This file contains the LocalConvection Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_LOCAL_CONVECTION_BOUNDARY_3D_HH
#define SET_LOCAL_CONVECTION_BOUNDARY_3D_HH

#include "setLocalConvectionBoundary3D.h"

namespace olb {

///Initialising the setLocalConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry<T,3>& superGeometry, int material,
                                T* uAv)
{
  setLocalConvectionBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry.getMaterialIndicator(material), uAv);
}

///Initialising the setLocalConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T* uAv)
{
  OstreamManager clout(std::cout, "setLocalConvectionBoundary");
  int _overlap = 0;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setLocalConvectionBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc),omega, indicator->getBlockIndicatorF(iCloc),
        uAv, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Add convection boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(BlockLattice<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, T* uAv, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

      if (discreteNormal[0] == 0) {//set postProcessors for indicated boundary cells
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          postProcessor = nullptr;
        }
        if (postProcessor) {
          _block.addPostProcessor(*postProcessor);
        }
      }
    }
  });
}

}//namespace olb
#endif

