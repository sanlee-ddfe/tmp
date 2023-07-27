/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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


#ifndef PARTICLE_PARALLEL_CREATOR_FUNCTIONS_3D_H
#define PARTICLE_PARALLEL_CREATOR_FUNCTIONS_3D_H

#include <algorithm>

#include "particles/functions/lambdaLoops.h"
#include "particles/functions/particleCreatorHelperFunctions.h"

namespace olb {

namespace particles {

namespace creators {

//// SUBGRID

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addSubgridSphere3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T radius, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.) )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidGeometry = sParticleSystem.getSuperStructure().getCuboidGeometry();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = cuboidGeometry.getC( position, globiCcentre );
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }

  //Iterate over particle systems
  std::size_t localID=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      addSubgridObject( particleSystem,
        position, radius, density, velocity );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Set localID
      localID=particle.getId();
    }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localID);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addSubgridSphereWithSpecies3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T radius, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const int species=0 )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidGeometry = sParticleSystem.getSuperStructure().getCuboidGeometry();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = cuboidGeometry.getC( position, globiCcentre );
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }

  //Iterate over particle systems
  std::size_t localID=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      addSubgridObject( particleSystem,
        position, radius, density, velocity );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Add species type to particle
      if constexpr( access::providesSpecies<PARTICLETYPE>() ){
        particle.template setField<PHYSPROPERTIES,SPECIES>( species );
      }
      //Set localID
      localID=particle.getId();
    }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localID);
}






} //namespace creators

} //namespace particles

} //namespace olb

#endif
