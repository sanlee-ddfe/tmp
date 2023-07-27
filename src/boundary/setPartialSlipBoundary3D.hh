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

//This file contains the Partial Slip Boundary
//This is an onLattice boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_PARTIAL_SLIP_BOUNDARY_HH
#define SET_PARTIAL_SLIP_BOUNDARY_HH

#include "setPartialSlipBoundary3D.h"

namespace olb {

///Initialising the PartialSlipBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary( SuperLattice<T, DESCRIPTOR>& sLattice, T tuner, SuperGeometry<T,3>& superGeometry, int material)
{
  setPartialSlipBoundary<T,DESCRIPTOR>(sLattice, tuner, superGeometry.getMaterialIndicator(material));
}

///Initialising the PartialSlipBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T tuner, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = 1;
  OstreamManager clout(std::cout, "setPartialslipBoundary");

  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setPartialSlipBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc),
                                         tuner, indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Set Partial Slip boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(BlockLattice<T,DESCRIPTOR>& block,T tuner, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setPartialSlipBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      if (tuner < 0. || tuner > 1.) {
        clout << "Warning: Could not setPartialSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), tuner must be between 0.1 and instead is=" << tuner <<", set to bounceBack" << std::endl;
        block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
      else {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {//set PostProcessor on indicated cells
          OstreamManager clout(std::cout, "setPartialslipBoundary");
          bool _output = false;
          if (_output) {
            clout << "setPartialSlipBoundary<" << discreteNormal[1] << ","<< discreteNormal[2] << ","<< discreteNormal[3] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << ", " << iZ << ", " << iZ << " )" << std::endl;
          }
          PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new PartialSlipBoundaryProcessorGenerator3D<T, DESCRIPTOR>(tuner, iX, iX, iY, iY, iZ, iZ, discreteNormal[1], discreteNormal[2], discreteNormal[3]);
          if (postProcessor) {
            block.addPostProcessor(*postProcessor);
          }
        }
        else {//define Dynamics on indicated cells
          clout << "Warning: Could not setPartialSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
          block.template defineDynamics<BounceBack>({iX, iY, iZ});
        }
      }
    }
  });
}


}//namespace olb
#endif
