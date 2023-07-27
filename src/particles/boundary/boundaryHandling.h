/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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



#ifndef PARTICLE_BOUNDARY_HANDLING_H
#define PARTICLE_BOUNDARY_HANDLING_H


namespace olb {

namespace particles {

namespace boundaries {


//Return normal on closest surface specified by solid boundary
//-referenceLength: For particles, use circumRadius etc.
template<typename T, unsigned D>
Vector<T,D> getNormalOnClosestSurface( SolidBoundary<T,D>& solidBoundary,
  Vector<T,D>& position, T referenceLength )
{
  T delXofCentralDifference = 1.e-8*referenceLength;
  return solidBoundary.getIndicator()->surfaceNormal(position, delXofCentralDifference);
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getNormalOnClosestSurface( SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  Particle<T,PARTICLETYPE>& particle )
{
  auto position = access::getPosition( particle );
  auto radius = access::getRadius( particle );
  return getNormalOnClosestSurface( solidBoundary, position, radius );
}

/// Generic treatment of particle wall contact
/// - evaluation of penetration by
///   - Version A: signed distance function
///   - Version B: simple min/max cubic bounds evaluation
/// - can be used for:
///   - f( particle, normal, penetrationDepth )
///   - f( particle, normal )
///   - f( particle )
template<bool useCubicBounds=false, typename T, typename PARTICLETYPE, typename F>
void doAtParticleWallContact(
  Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  F f )
{
  using namespace descriptors;
  //Retrieve quantities
  auto position = access::getPosition( particle );
  T radius = access::getRadius( particle );
  //VERSION A: Use signed distance functions of indicator
  if constexpr (!useCubicBounds){
    T distToWall = solidBoundary.getIndicator()->signedDistance(position)-radius;
    if (distToWall<=0){

      if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&,T>){
        auto normal = getNormalOnClosestSurface( solidBoundary, position, radius );
        T penetrationDepth = -distToWall;
        f( particle, normal, penetrationDepth );
      } else if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&>){
        auto normal = getNormalOnClosestSurface( solidBoundary, position, radius );
        f( particle, normal );
      } else {
        f( particle );
      }

    }
  }
  //VERSION B: Use origin and end of cubic indicator bounds (min, max)
  else {
    auto origin = solidBoundary.getIndicator()->getMin();
    auto end = solidBoundary.getIndicator()->getMax();
    for (int iDim=0; iDim<PARTICLETYPE::d; ++iDim){
      T posParticleStart = position[iDim] - radius;
      T posParticleEnd = position[iDim] + radius;
      //Retrieve distance for lower and upper limit
      T distToWallOrig = posParticleStart - origin[iDim];
      T distToWallEnd = end[iDim] - posParticleEnd;
      //Set up normal
      Vector<T,PARTICLETYPE::d> normal(0.);
      //Evaluate distances
      if (distToWallOrig<=0){
        normal[iDim] = 1;
        if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&,T>){
          T penetrationDepth = -distToWallOrig;
          f( particle, normal, penetrationDepth );
        } else if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&>){
          f( particle, normal );
        } else {
          f( particle );
        }
      } else if (distToWallEnd<=0){
        normal[iDim] = -1;
        if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&,T>){
          T penetrationDepth = -distToWallEnd;
          f( particle, normal, penetrationDepth );
        } else if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&,Vector<T,PARTICLETYPE::d>&>){
          f( particle, normal );
        } else {
          f( particle );
        }
      }
    }
  } //constexpr if (!useCubicBounds)
}


//Orientation between given normal and closest surface of solid boundary
// - usable for quantifying tilting of surface attached particles
template<typename T, unsigned D>
Vector<T,utilities::dimensions::convert<D>::rotation> getRelativeSurfaceOrientation(
  SolidBoundary<T,D>& solidBoundary, Vector<T,D>& position,
  Vector<T,D>& normalToBeComparedTo, T referenceLength )
{
  auto normalWall = getNormalOnClosestSurface( solidBoundary, position, referenceLength );
  auto orientation = util::angleBetweenVectors( normalWall, normalToBeComparedTo );
  return orientation;
}

} //namespace boundaries

} //namespace particles

} //namespace olb


#endif
