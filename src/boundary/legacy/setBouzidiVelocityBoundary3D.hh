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

//This file contains the BouzidiVelocityBoundary
//This is an offLattice Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_BOUZIDI_VELOCITY_BOUNDARY_HH
#define SET_BOUZIDI_VELOCITY_BOUNDARY_HH

#include "setBouzidiVelocityBoundary3D.h"

namespace olb {

namespace legacy {

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the BouzidiVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,SuperGeometry<T,3>& superGeometry, int material,
                                IndicatorF3D<T>& geometryIndicator,
                                std::vector<int> bulkMaterials)
{

  setBouzidiVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
      superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
      geometryIndicator);
}

///Initialising the BouzidiVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                                FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
                                IndicatorF3D<T>&                   geometryIndicator)
{
  OstreamManager clout(std::cout, "setBouzidiBoundary");
  T _epsFraction = 0.0001;
  int overlap = 1;
  bool _output = false;
  if (_output) {
    clout << "epsFraction=" << _epsFraction << std::endl;
    clout.setMultiOutput(true);
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    if (_output) {
      clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
            << " starts to read distances for Velocity Boundary..." << std::endl;
    }
    //this is a guess to solve the private member block issue
    setBouzidiVelocityBoundary<T,DESCRIPTOR>(
      sLattice.getBlock(iCloc), boundaryIndicator->getBlockIndicatorF(iCloc),
      bulkIndicator->getBlockIndicatorF(iCloc),
      geometryIndicator, _epsFraction);
    if (_output) {
      clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
            << " finished reading distances for Velocity Boundary." << std::endl;
    }
  }
  if (_output) {
    clout.setMultiOutput(false);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice,std::forward<decltype(boundaryIndicator)>(boundaryIndicator), overlap);
}


////////// BlockLattice Domain  /////////////////////////////////////////


template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockIndicatorF3D<T>& boundaryIndicator,
                                BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator, T _epsFraction)
{
  block.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (boundaryIndicator(iX,iY,iZ)) {
      setBouzidiVelocityBoundary<T,DESCRIPTOR>(block,
                                               bulkIndicator.getBlockGeometry(),
                                               iX, iY, iZ,
                                               geometryIndicator,
                                               bulkIndicator,
                                               _epsFraction);
    }
  });
}


template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int iX, int iY, int iZ,
                                IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator, T _epsFraction)
{
  OstreamManager clout(std::cout, "setBouzidiBoundary");
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    const int iZn = iZ + c[2];
    if (blockGeometryStructure.isInside(iXn,iYn,iZn) && bulkIndicator(iXn,iYn,iZn)) {
      T dist = -1;
      T physR[3];
      blockGeometryStructure.getPhysR(physR,{iXn,iYn,iZn});
      T voxelSize=blockGeometryStructure.getDeltaR();

      Vector<T,3> physC(physR);
      Vector<T,3> direction(-voxelSize*c[0],-voxelSize*c[1],-voxelSize*c[2]);
      T cPhysNorm = voxelSize*util::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

      if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
        T epsX = voxelSize*c[0]*_epsFraction;
        T epsY = voxelSize*c[1]*_epsFraction;
        T epsZ = voxelSize*c[2]*_epsFraction;

        Vector<T,3> physC2(physC);
        physC2[0] += epsX;
        physC2[1] += epsY;
        physC2[2] += epsZ;
        Vector<T,3> direction2(direction);
        direction2[0] -= 2.*epsX;
        direction2[1] -= 2.*epsY;
        direction2[2] -= 2.*epsZ;

        if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
          clout << "ERROR: no boundary found at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << descriptors::opposite<DESCRIPTOR >(iPop)
                << std::endl;

        }
        T distNew = (dist - util::sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
        if (distNew < 0.5) {
          dist = 0;
        }
        else {
          dist = 0.5 * cPhysNorm;
          clout << "WARNING: distance at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << descriptors::opposite<DESCRIPTOR >(iPop) << ": "
                << distNew
                << " rounded to "
                << dist/cPhysNorm
                << std::endl;

        }
      }
      distances[descriptors::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
    } // bulk indicator if
  } // iPop loop
  setBouzidiVelocityBoundary<T,DESCRIPTOR>(block, blockGeometryStructure, iX, iY, iZ, distances);
}


template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q])
{
  T location[DESCRIPTOR::d];
  blockGeometryStructure.getPhysR(location, {x,y,z});
  T distancesCopy[DESCRIPTOR::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  Dynamics<T,DESCRIPTOR>* dynamics = new legacy::OffDynamics<T, DESCRIPTOR>(location, distancesCopy);
  block.defineDynamics({x,y,z}, dynamics);

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    if (!util::nearZero(distances[iPop]+1)) {
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);

      setBouzidiVelocityBoundary<T,DESCRIPTOR>(block, blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop]);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void setBouzidiVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int x, int y, int z, int iPop, T dist)
{
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  //setOnePointVelocityBoundary. This corresponds to Bounce Back
  auto postProcessor = std::unique_ptr<PostProcessorGenerator3D<T, DESCRIPTOR>> { nullptr };
  if (blockGeometryStructure.getMaterial({x-c[0], y-c[1], z-c[2]}) != 1) {
    postProcessor.reset(new VelocityBounceBackPostProcessorGenerator3D <T, DESCRIPTOR>(x, y, z, iPop, dist));
  }
  //setTwoPointVelocityBoundary.This corresponds to Linear Interpolation
  else {
    postProcessor.reset(new VelocityBouzidiLinearPostProcessorGenerator3D<T, DESCRIPTOR>(x, y, z, iPop, dist));
  }
  if (postProcessor && !block.isPadding({x,y,z})) {
    block.addPostProcessor(*postProcessor);
  }
}

}

}

#endif
