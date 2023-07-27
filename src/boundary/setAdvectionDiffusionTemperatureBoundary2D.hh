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

//This file contains the AdvectionDiffusionTemperature Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_ADVECTION_DIFFUSION_TEMPERATURE_BOUNDARY_2D_HH
#define SET_ADVECTION_DIFFUSION_TEMPERATURE_BOUNDARY_2D_HH

#include "setAdvectionDiffusionTemperatureBoundary2D.h"


namespace olb {

///Initialising the AdvectionDiffusionTemperatureBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega,SuperGeometry<T,2>& superGeometry, int material, bool neumann)
{
  setAdvectionDiffusionTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material), neumann);
}

///Initialising the AdvectionDiffusionTemperatureBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega,FunctorPtr<SuperIndicatorF2D<T>>&& indicator, bool neumann)
{
  OstreamManager clout(std::cout, "setAdvectionDiffusionTemperatureBoundary");

  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setAdvectionDiffusionTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc), omega,
        indicator->getBlockIndicatorF(iCloc), includeOuterCells, neumann);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);

}

///set AdvectionDiffusionTemperature boundary on indicated cells inside the block domains
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setAdvectionDiffusionTemperatureBoundary(BlockLattice<T,DESCRIPTOR>& block, T omega, BlockIndicatorF2D<T>& indicator,
    bool includeOuterCells, bool neumann)
{
  using namespace boundaryhelper;
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(3,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX,iY);
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;

      if (neumann){
        block.addPostProcessor(
          typeid(stage::PostStream), {iX,iY},
          promisePostProcessorForNormal<T, DESCRIPTOR, NeumannBoundarySingleLatticePostProcessor2D>(
            Vector <int,2> (discreteNormal.data()+1)
          )
        );
      }

      if (discreteNormal[0] == 0) {
          dynamics = block.getDynamics(DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
            AdvectionDiffusionBoundariesDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
      }
      else if (discreteNormal[0] == 1) {//set momenta, dynamics and postProcessors for indicated AdvectionDiffusionTemperatureBoundary corner cells
          dynamics = block.getDynamics(NormalMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
            AdvectionDiffusionCornerDynamics2D,MixinDynamics,momenta::EquilibriumBoundaryTuple
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
      }
      dynamics->getParameters(block).template set<descriptors::OMEGA>(omega);
      setBoundary(block, iX,iY, dynamics);
    }
  });
}
}//namespace olb


#endif
