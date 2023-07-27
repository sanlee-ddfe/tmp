/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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

/** \file
 * The description of a vector of 3D cuboid  -- generic implementation.
 */


#ifndef CUBOID_GEOMETRY_3D_HH
#define CUBOID_GEOMETRY_3D_HH


#include <iostream>
#include <math.h>
#include <algorithm>
#include <set>
#include <limits>

#include "geometry/cuboidGeometry3D.h"
#include "geometry/cuboidGeometryMinimizer.h"

namespace olb {

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D()
  : _motherCuboid(0,0,0,0,0,0,0), _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  _cuboids.reserve(2);
  add(_motherCuboid);
  split(0, 1);
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(T originX, T originY, T originZ, T deltaR,
                                      int nX, int nY, int nZ, int nC)
  : _motherCuboid(originX, originY, originZ, deltaR, nX, nY, nZ),
    _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  _cuboids.reserve(nC+2);
  add(_motherCuboid);
  split(0, nC);
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(std::vector<T> origin, T deltaR,
                                      std::vector<int> extent, int nC)
  : CuboidGeometry3D(origin[0], origin[1], origin[2], deltaR,
                     extent[0], extent[1], extent[2], nC)
{ }

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(Vector<T,3> origin, T deltaR,
                                      Vector<int,3> extent, int nC)
  : CuboidGeometry3D(origin[0], origin[1], origin[2], deltaR,
                     extent[0], extent[1], extent[2], nC)
{ }

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(const Cuboid3D<T>& motherCuboid, int nC)
  : CuboidGeometry3D(
    motherCuboid.getOrigin(), motherCuboid.getDeltaR(),
    motherCuboid.getExtent(), nC)
{ }

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(IndicatorF3D<T>& indicatorF, T voxelSize, int nC)
  : _motherCuboid( indicatorF.getMin()[0],  indicatorF.getMin()[1], indicatorF.getMin()[2], voxelSize,
                   (int)((indicatorF.getMax()[0] - indicatorF.getMin()[0]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[1] - indicatorF.getMin()[1]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[2] - indicatorF.getMin()[2]) / voxelSize + 1.5)),
    _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  _cuboids.reserve(nC+2);
  add(_motherCuboid);
  if (nC > 1) {
    split(0, nC);
    shrink(indicatorF);
  }
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(std::shared_ptr<IndicatorF3D<T>> indicator_sharedPtrF, T voxelSize, int nC)
  : CuboidGeometry3D<T>(*indicator_sharedPtrF, voxelSize, nC)
{
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(IndicatorF3D<T>& indicatorF, T voxelSize, int nC, std::string minimizeBy)
  : _motherCuboid( indicatorF.getMin()[0],  indicatorF.getMin()[1], indicatorF.getMin()[2], voxelSize,
                   (int)((indicatorF.getMax()[0] - indicatorF.getMin()[0]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[1] - indicatorF.getMin()[1]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[2] - indicatorF.getMin()[2]) / voxelSize + 1.5)),
    _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  _cuboids.reserve(nC+2);
  add(_motherCuboid);

  if ( minimizeBy == "volume" ) {
    minimizeByVolume(*this, indicatorF, nC);
  }
  else if ( minimizeBy == "weight" ) {
    minimizeByWeight(*this, indicatorF, nC);
  }
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(std::shared_ptr<IndicatorF3D<T>> indicator_sharedPtrF, T voxelSize, int nC, std::string minimizeBy)
  : CuboidGeometry3D<T>(*indicator_sharedPtrF, voxelSize, nC, minimizeBy)
{
}

template<typename T>
CuboidGeometry3D<T>::~CuboidGeometry3D() {};


template<typename T>
Cuboid3D<T>& CuboidGeometry3D<T>::get(int iC)
{
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> const& CuboidGeometry3D<T>::get(int iC) const
{
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> CuboidGeometry3D<T>::getMotherCuboid()
{
  return _motherCuboid;
}

template<typename T>
Cuboid3D<T> const& CuboidGeometry3D<T>::getMotherCuboid() const
{
  return _motherCuboid;
}

template<typename T>
void CuboidGeometry3D<T>::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
  _periodicityOn[0] = periodicityX;
  _periodicityOn[1] = periodicityY;
  _periodicityOn[2] = periodicityZ;
}


template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z, int offset) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z, offset)) {
      return (int)i;
    }
  }
  return (int)i;
}


template<typename T>
int CuboidGeometry3D<T>::get_iC(Vector<T,3> coords, int offset) const
{
  return get_iC(coords[0], coords[1], coords[2], offset);
}

template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z, int orientationX, int orientationY,
                                int orientationZ) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z) &&
        _cuboids[i].checkPoint(x + orientationX / _cuboids[i].getDeltaR(),
                               y + orientationY / _cuboids[i].getDeltaR(),
                               z + orientationZ / _cuboids[i].getDeltaR())) {
      return (int)i;
    }
  }
  return (int)i;
}

template<typename T>
bool CuboidGeometry3D<T>::getC(T physR[3], int& iC) const
{
  const int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getC(std::vector<T> physR, int& iC) const
{
  const int iCtmp = get_iC(physR);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getC(const Vector<T,3>& physR, int& iC) const
{
  const int iCtmp = get_iC(physR);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getLatticeR(int latticeR[4], const T physR[3]) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)util::floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[2] = (int)util::floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[3] = (int)util::floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)util::floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)util::floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[3] = (int)util::floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] ) / _cuboids[latticeR[0]].getDeltaR() );
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getFloorLatticeR(
  const Vector<T,3>& physR, Vector<int,4>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)util::floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)util::floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[3] = (int)util::floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] ) / _cuboids[latticeR[0]].getDeltaR() );
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], const int& iCglob, const int& iX, const int& iY, const int& iZ) const
{
  _cuboids[iCglob].getPhysR(physR, iX, iY, iZ);
  for (int iDim = 0; iDim < 3; iDim++) {
    if (_periodicityOn[iDim]) {
      //std::cout << "!!! " << iDim << _periodicityOn[iDim] <<":"<< _motherCuboid.getDeltaR()*(_motherCuboid.getExtent()[iDim]) << std::endl;
      physR[iDim] = remainder( physR[iDim] - _motherCuboid.getOrigin()[iDim]
                               + _motherCuboid.getDeltaR() * (_motherCuboid.getExtent()[iDim]),
                               _motherCuboid.getDeltaR() * (_motherCuboid.getExtent()[iDim]));
      // solving the rounding error problem for double
      if ( physR[iDim]*physR[iDim] < 0.001 * _motherCuboid.getDeltaR()*_motherCuboid.getDeltaR() ) {
        physR[iDim] = T();
      }
      // make it to mod instead remainer
      if ( physR[iDim] < 0 ) {
        physR[iDim] += _motherCuboid.getDeltaR() *( _motherCuboid.getExtent()[iDim]);
      }
      // add origin
      physR[iDim] += _motherCuboid.getOrigin()[iDim];
    }
  }
  return;
}

template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], const int latticeR[4]) const
{
  getPhysR(physR, latticeR[0],  latticeR[1],  latticeR[2],  latticeR[3]);
}

template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], LatticeR<4> latticeR) const
{
  getPhysR(physR, latticeR[0],  latticeR[1],  latticeR[2],  latticeR[3]);
}

template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], const int iCglob, LatticeR<3> latticeR) const
{
  getPhysR(physR, iCglob,  latticeR[0],  latticeR[1],  latticeR[2]);
}

template<typename T>
void CuboidGeometry3D<T>::getNeighbourhood(int cuboid, std::vector<int>& neighbours, int overlap)
{
  neighbours.clear();

  std::set<int> dummy;

  for (int iC = 0; iC < getNc(); iC++) {
    if (cuboid == iC) {
      continue;
    }
    T globX = get(iC).getOrigin()[0];
    T globY = get(iC).getOrigin()[1];
    T globZ = get(iC).getOrigin()[2];
    T nX = get(iC).getNx();
    T nY = get(iC).getNy();
    T nZ = get(iC).getNz();
    T deltaR = get(iC).getDeltaR();
    if (get(cuboid).checkInters(globX,
                                globX + (nX + overlap - 1)*deltaR,
                                globY - overlap * deltaR,
                                globY + (nY + overlap - 1)*deltaR,
                                globZ - overlap * deltaR,
                                globZ + (nZ + overlap - 1)*deltaR, overlap)) {
      //neighbours.push_back(iC);
      dummy.insert(iC);
    }

    if (_periodicityOn[0]) {
      if (get(cuboid).getOrigin()[0] + (get(cuboid).getNx() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[0]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0]-getMaxPhysR()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[0] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[0]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0]+getMaxPhysR()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

    if (_periodicityOn[1]) {
      if (get(cuboid).getOrigin()[1] + (get(cuboid).getNy() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[1]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1]-getMaxPhysR()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[1] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[1]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1]+getMaxPhysR()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

    if (_periodicityOn[2]) {
      if (get(cuboid).getOrigin()[2] + (get(cuboid).getNz() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[2]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2]-getMaxPhysR()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[2] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[2]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2]+getMaxPhysR()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

  }
  std::set<int>::iterator it = dummy.begin();
  for (; it != dummy.end(); ++it) {
    neighbours.push_back(*it);
  }
}

template<typename T>
int CuboidGeometry3D<T>::getNc() const
{
  return _cuboids.size();
}

template<typename T>
T CuboidGeometry3D<T>::getMinRatio() const
{
  T minRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() < minRatio) {
      minRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() < minRatio) {
      minRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
    }
    if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() < minRatio) {
      minRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
    }
  }
  return minRatio;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxRatio() const
{
  T maxRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
    }
    if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
    }
  }
  return maxRatio;
}

template<typename T>
Vector<T,3> CuboidGeometry3D<T>::getMinPhysR() const
{
  Vector<T,3> output (_cuboids[0].getOrigin());
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] < output[0]) {
      output[0] = _cuboids[i].getOrigin()[0];
    }
    if (_cuboids[i].getOrigin()[1] < output[1]) {
      output[1] = _cuboids[i].getOrigin()[1];
    }
    if (_cuboids[i].getOrigin()[2] < output[2]) {
      output[2] = _cuboids[i].getOrigin()[2];
    }
  }
  return output;
}

template<typename T>
Vector<T,3> CuboidGeometry3D<T>::getMaxPhysR() const
{
  Vector<T,3> output (_cuboids[0].getOrigin());
  output[0] += _cuboids[0].getNx()*_cuboids[0].getDeltaR();
  output[1] += _cuboids[0].getNy()*_cuboids[0].getDeltaR();
  output[2] += _cuboids[0].getNz()*_cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR() > output[0]) {
      output[0] = _cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR();
    }
    if (_cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR() > output[1]) {
      output[1] = _cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR();
    }
    if (_cuboids[i].getOrigin()[2] + _cuboids[i].getNz()*_cuboids[i].getDeltaR() > output[2]) {
      output[2] = _cuboids[i].getOrigin()[2] + _cuboids[i].getNz()*_cuboids[i].getDeltaR();
    }
  }
  return output;
}

template<typename T>
T CuboidGeometry3D<T>::getMinPhysVolume() const
{
  T minVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() < minVolume) {
      minVolume = _cuboids[i].getPhysVolume();
    }
  }
  return minVolume;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxPhysVolume() const
{
  T maxVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() > maxVolume) {
      maxVolume = _cuboids[i].getPhysVolume();
    }
  }
  return maxVolume;
}

template<typename T>
size_t CuboidGeometry3D<T>::getMinLatticeVolume() const
{
  size_t minNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() < minNodes) {
      minNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return minNodes;
}

template<typename T>
size_t CuboidGeometry3D<T>::getMaxLatticeVolume() const
{
  size_t maxNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() > maxNodes) {
      maxNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return maxNodes;
}

template<typename T>
size_t CuboidGeometry3D<T>::getMinLatticeWeight() const
{
  size_t minNodes = _cuboids[0].getWeight();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getWeight() < minNodes) {
      minNodes = _cuboids[i].getWeight();
    }
  }
  return minNodes;
}

template<typename T>
size_t CuboidGeometry3D<T>::getMaxLatticeWeight() const
{
  size_t maxNodes = _cuboids[0].getWeight();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getWeight() > maxNodes) {
      maxNodes = _cuboids[i].getWeight();
    }
  }
  return maxNodes;
}

template<typename T>
T CuboidGeometry3D<T>::getMinDeltaR() const
{
  T minDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() < minDelta) {
      minDelta = _cuboids[i].getDeltaR();
    }
  }
  return minDelta;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxDeltaR() const
{
  T maxDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() > maxDelta) {
      maxDelta = _cuboids[i].getDeltaR();
    }
  }
  return maxDelta;
}


template<typename T>
bool CuboidGeometry3D<T>::operator==(CuboidGeometry3D<T>& rhs)
{
  return     _motherCuboid == rhs._motherCuboid
             && _periodicityOn == rhs._periodicityOn
             &&       _cuboids == rhs._cuboids;
}



template<typename T>
void CuboidGeometry3D<T>::add(Cuboid3D<T> cuboid)
{

  _cuboids.push_back(cuboid);
}

template<typename T>
void CuboidGeometry3D<T>::remove(int iC)
{

  _cuboids.erase(_cuboids.begin() + iC);
}


template<typename T>
void CuboidGeometry3D<T>::remove(IndicatorF3D<T>& indicatorF)
{
  std::vector<bool> allZero;
  int latticeR[4];
  T physR[3];
  for (unsigned iC = 0; iC < _cuboids.size(); iC++) {
    latticeR[0] = iC;
    allZero.push_back(true);
    for (int iX = 0; iX < _cuboids[iC].getNx(); iX++) {
      for (int iY = 0; iY < _cuboids[iC].getNy(); iY++) {
        for (int iZ = 0; iZ < _cuboids[iC].getNz(); iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          getPhysR(physR,latticeR);
          bool inside[1];
          indicatorF(inside,physR);
          if (inside[0]) {
            allZero[iC] = 0;
          }
        }
      }
    }
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (allZero[iC] ) {
      remove(iC);
    }
  }
}

template<typename T>
void CuboidGeometry3D<T>::removeByWeight()
{
  std::vector<bool> allZero(_cuboids.size(), false);
  for (unsigned iC = 0; iC < _cuboids.size(); iC++) {
    allZero[iC] = (_cuboids[iC].getWeight() == 0);
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (allZero[iC]) {
      remove(iC);
    }
  }
}

template<typename T>
void CuboidGeometry3D<T>::shrink(int iC, IndicatorF3D<T>& indicatorF)
{
  int latticeR[4];
  T physR[3];
  bool inside[1];

  latticeR[0] = iC;
  size_t fullCells = 0;
  int xN = get(iC).getNx();
  int yN = get(iC).getNy();
  int zN = get(iC).getNz();
  int maxX = 0;
  int maxY = 0;
  int maxZ = 0;
  int newX = xN - 1;
  int newY = yN - 1;
  int newZ = zN - 1;
  for (int iX = 0; iX < xN; iX++) {
    for (int iY = 0; iY < yN; iY++) {
      for (int iZ = 0; iZ < zN; iZ++) {
        latticeR[1] = iX;
        latticeR[2] = iY;
        latticeR[3] = iZ;
        getPhysR(physR,latticeR);
        indicatorF(inside,physR);
        if (inside[0]) {
          fullCells++;
          maxX = util::max(maxX, iX);
          maxY = util::max(maxY, iY);
          maxZ = util::max(maxZ, iZ);
          newX = util::min(newX, iX);
          newY = util::min(newY, iY);
          newZ = util::min(newZ, iZ);
        }
      }
    }
  }
  //    if (maxX+2 < xN) maxX+=2; else if (maxX+1 < xN) maxX+=1;
  //    if (maxY+2 < yN) maxY+=2; else if (maxY+1 < yN) maxY+=1;
  //    if (maxZ+2 < zN) maxZ+=2; else if (maxZ+1 < zN) maxZ+=1;
  //
  //    if (newX-2 >= 0) newX-=2; else if (newX-1 >= 0) newX-=1;
  //    if (newY-2 >= 0) newY-=2; else if (newY-1 >= 0) newY-=1;
  //    if (newZ-2 >= 0) newZ-=2; else if (newZ-1 >= 0) newZ-=1;

  if (fullCells > 0) {
    get(iC).setWeight(fullCells);
    _cuboids[iC].resize(newX, newY, newZ, maxX - newX + 1, maxY - newY + 1, maxZ - newZ + 1);
  }
  else {
    remove(iC);
  }
}

template<typename T>
void CuboidGeometry3D<T>::shrink(IndicatorF3D<T>& indicatorF)
{
  for (int iC = getNc() - 1; iC >= 0; iC--) {
    shrink(iC, indicatorF);
  }
  // shrink mother cuboid
  Vector<T,3> minPhysR = getMinPhysR();
  Vector<T,3> maxPhysR = getMaxPhysR();
  T minDelataR = getMinDeltaR();
  _motherCuboid = Cuboid3D<T>(minPhysR[0], minPhysR[1], minPhysR[2], minDelataR,
                              (int)((maxPhysR[0]-minPhysR[0])/minDelataR + 0.5),
                              (int)((maxPhysR[1]-minPhysR[1])/minDelataR + 0.5),
                              (int)((maxPhysR[2]-minPhysR[2])/minDelataR + 0.5));
}


template<typename T>
void CuboidGeometry3D<T>::split(int iC, int p)
{
  Cuboid3D<T> temp(_cuboids[iC].getOrigin()[0], _cuboids[iC].getOrigin()[1],
                   _cuboids[iC].getOrigin()[2],  _cuboids[iC].getDeltaR(),
                   _cuboids[iC].getNx(), _cuboids[iC].getNy(), _cuboids[iC].getNz());
  temp.divide(p, _cuboids);
  remove(iC);
}

template<typename T>
void CuboidGeometry3D<T>::splitRegular(int iC, int width)
{
  Cuboid3D<T> temp(_cuboids[iC].getOrigin()[0], _cuboids[iC].getOrigin()[1],
                   _cuboids[iC].getOrigin()[2],  _cuboids[iC].getDeltaR(),
                   _cuboids[iC].getNx(), _cuboids[iC].getNy(), _cuboids[iC].getNz());
  const int p = std::max(1, temp.getNx() / width);
  const int q = std::max(1, temp.getNy() / width);
  const int r = std::max(1, temp.getNz() / width);
  temp.divide(p, q, r, _cuboids);
  remove(iC);
}

template<typename T>
void CuboidGeometry3D<T>::splitByWeight(int iC, int p, IndicatorF3D<T>& indicatorF)
{
  T averageWeight = get(iC).getWeight() / (T) p;
  // clout << "Mother " << get(iC).getWeight() << " " << averageWeight << std::endl;
  Cuboid3D<T> temp(_cuboids[iC].getOrigin()[0], _cuboids[iC].getOrigin()[1],
                   _cuboids[iC].getOrigin()[2],  _cuboids[iC].getDeltaR(),
                   _cuboids[iC].getNx(), _cuboids[iC].getNy(), _cuboids[iC].getNz());

  int latticeR[4];
  T physR[3];
  latticeR[0] = iC;
  int xN = get(iC).getNx();
  int yN = get(iC).getNy();
  int zN = get(iC).getNz();
  T deltaR = get(iC).getDeltaR();
  int fullCells = 0;

  Vector<T, 3> globPos_child = get(iC).getOrigin();
  std::vector<int> extend_child = {xN, yN, zN};
  int localPos_child = 0;

  // looking for largest extend, because halfing the cuboid by its largest extend will result in the smallest surface and therfore in the least comminication cells
  if ( get(iC).getNx() >= get(iC).getNy() && get(iC).getNx() >= get(iC).getNz()) {
    // clout << "Cut in x direction!" << std::endl;

    // for each child cuboid search for the optimal cutting plane
    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iX = localPos_child; iX < xN; iX++) {
        fullCells_minusOne = fullCells;
        for (int iY = 0; iY < yN; iY++) {
          for (int iZ = 0; iZ < zN; iZ++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            getPhysR(physR,latticeR);
            bool inside[1];
            indicatorF(inside,physR);
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        // the optimal cutting plane is determined, so that the child cuboid's cells inside the indicator are the closest to the total cells inside the indicator per number of children
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iX--;
          }
          // clout << "found optimal iX = " << iX << std::endl;
          extend_child[0] = iX - localPos_child + 1;

          Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
          _cuboids.push_back(child);

          globPos_child[0] += extend_child[0]*deltaR;
          localPos_child += extend_child[0] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[0] = xN - localPos_child + p - 1;

    Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;

  }
  else if ( get(iC).getNy() >= get(iC).getNx() && get(iC).getNy() >= get(iC).getNz()) {
    // clout << "Cut in y direction!" << std::endl;

    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iY = localPos_child; iY < yN; iY++) {
        fullCells_minusOne = fullCells;
        for (int iX = 0; iX < xN; iX++) {
          for (int iZ = 0; iZ < zN; iZ++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            getPhysR(physR,latticeR);
            bool inside[1];
            indicatorF(inside,physR);
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iY--;
          }
          // clout << "found optimal iY = " << iY << std::endl;
          extend_child[1] = iY - localPos_child + 1;

          Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
          _cuboids.push_back(child);

          globPos_child[1] += extend_child[1]*deltaR;
          localPos_child += extend_child[1] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[1] = yN - localPos_child + p - 1;

    Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;
  }
  else {
    // clout << "Cut in z direction!" << std::endl;

    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iZ = localPos_child; iZ < zN; iZ++) {
        fullCells_minusOne = fullCells;
        for (int iY = 0; iY < yN; iY++) {
          for (int iX = 0; iX < xN; iX++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            getPhysR(physR,latticeR);
            bool inside[1];
            indicatorF(inside,physR);
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iZ--;
          }
          // clout << "found optimal iZ = " << iZ << std::endl;
          extend_child[2] = iZ - localPos_child + 1;

          Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
          _cuboids.push_back(child);

          globPos_child[2] += extend_child[2]*deltaR;
          localPos_child += extend_child[2] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[2] = zN - localPos_child + p - 1;

    Cuboid3D<T> child(globPos_child[0], globPos_child[1], globPos_child[2], deltaR, extend_child[0], extend_child[1], extend_child[2]);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;
  }
}

template<typename T>
void CuboidGeometry3D<T>::splitFractional(int iC, int iD, std::vector<T> fractions)
{
  Cuboid3D<T> tmp = _cuboids[iC];
  tmp.divideFractional(iD, fractions, _cuboids);
  remove(iC);
}

template<typename T>
void CuboidGeometry3D<T>::swap(CuboidGeometry3D<T>& rhs)
{
  std::swap(this->_cuboids, rhs._cuboids);
  std::swap(this->_motherCuboid, rhs._motherCuboid);
  std::swap(this->_periodicityOn[0], rhs._periodicityOn[0]);
  std::swap(this->_periodicityOn[1], rhs._periodicityOn[1]);
  std::swap(this->_periodicityOn[2], rhs._periodicityOn[2]);
  std::swap(this->clout, rhs.clout);
}

template<typename T>
void CuboidGeometry3D<T>::swapCuboids(std::vector< Cuboid3D<T> >& cuboids)
{
  _cuboids.swap(cuboids);
}

template<typename T>
void CuboidGeometry3D<T>::replaceCuboids(std::vector< Cuboid3D<T> >& cuboids)
{
  this->_cuboids.clear();
  for ( unsigned iC = 0; iC < cuboids.size(); iC++) {
    add(cuboids[iC]);
  }
}

template<typename T>
std::size_t CuboidGeometry3D<T>::replaceContainedBy(Cuboid3D<T> mother)
{
  std::vector<int> toBeRemoved;
  for (unsigned iC=0; iC < _cuboids.size(); ++iC) {
    if (mother.partialOverlapWith(_cuboids[iC]) && !mother.contains(_cuboids[iC])) {
      return 0;
    } else if (mother.contains(_cuboids[iC])) {
      toBeRemoved.emplace_back(iC);
    }
  }
  std::sort(toBeRemoved.begin(), toBeRemoved.end(), std::greater{});
  mother.setWeight(0);
  for (int dC : toBeRemoved) {
    mother.setWeight(mother.getWeight() + _cuboids[dC].getWeight());
    _cuboids.erase(_cuboids.begin() + dC);
  }
  _cuboids.push_back(mother);
  return toBeRemoved.size();
}

template<typename T>
std::size_t CuboidGeometry3D<T>::hypotheticalReplaceContainedBy(Cuboid3D<T> mother)
{
  std::size_t weight = 0;
  for (unsigned iC=0; iC < _cuboids.size(); ++iC) {
    if (mother.partialOverlapWith(_cuboids[iC]) && !mother.contains(_cuboids[iC])) {
      return 0;
    } else if (mother.contains(_cuboids[iC])) {
      weight += _cuboids[iC].getWeight();
    }
  }
  return weight;
}

template<typename T>
void CuboidGeometry3D<T>::setWeights(IndicatorF3D<T>& indicatorF)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC=0; iC < getNc(); ++iC) {
    int latticeR[4] { iC, 0, 0, 0 };
    T physR[3];
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int zN = get(iC).getNz();
    size_t fullCells = 0;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        for (int iZ = 0; iZ < zN; iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          getPhysR(physR,latticeR);
          bool inside[1];
          indicatorF(inside,physR);
          if (inside[0]) {
            fullCells++;
          }
        }
      }
    }
    get(iC).setWeight(fullCells);
  }
}


template<typename T>
size_t CuboidGeometry3D<T>::getNblock() const
{
  return   1  // _periodicityOn
           + _motherCuboid.getNblock() // _motherCuboid
           + _cuboids.size() > 0 ? 1 + _cuboids.size() * _cuboids[0].getNblock() : 0; // _cuboids;
}


template<typename T>
size_t CuboidGeometry3D<T>::getSerializableSize() const
{
  return   3 * sizeof(bool)  // _periodicityOn
           + _motherCuboid.getSerializableSize() // _motherCuboid
           + (_cuboids.size() > 0 ?
              sizeof(size_t) + _cuboids.size() * _cuboids[0].getSerializableSize() :
              0); // _cuboids;
}

template<typename T>
bool* CuboidGeometry3D<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  size_t sizeBufferIndex = 0;
  bool* dataPtr = nullptr;

  registerVar<bool>                           (iBlock, sizeBlock, currentBlock, dataPtr, _periodicityOn[0], 3);
  registerSerializableOfConstSize             (iBlock, sizeBlock, currentBlock, dataPtr, _motherCuboid, loadingMode);
  registerStdVectorOfSerializablesOfConstSize (iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr,
      _cuboids, loadingMode);

  return dataPtr;
}


template<typename T>
void CuboidGeometry3D<T>::print() const
{
  clout << "---Cuboid Stucture Statistics---" << std::endl;
  clout << " Number of Cuboids: " << "\t" << getNc() << std::endl;
  clout << " Delta  (min): " << "\t" << "\t" << getMinDeltaR() << std::endl;
  clout << "        (max): " << "\t" << "\t" << getMaxDeltaR() << std::endl;
  clout << " Ratio  (min): " << "\t" << "\t" << getMinRatio() << std::endl;
  clout << "        (max): " << "\t" << "\t" << getMaxRatio() << std::endl;
  clout << " Nodes  (min): " << "\t" << "\t" << getMinLatticeVolume() << std::endl;
  clout << "        (max): " << "\t" << "\t" << getMaxLatticeVolume() << std::endl;
  clout << " Weight (min): " << "\t" << "\t" << getMinLatticeWeight() << std::endl;
  clout << "        (max): " << "\t" << "\t" << getMaxLatticeWeight() << std::endl;
  clout << "--------------------------------" << std::endl;
}

template<typename T>
void CuboidGeometry3D<T>::printExtended()
{
  clout << "Mothercuboid :" << std::endl;
  getMotherCuboid().print();

  for (int iC = 0; iC < getNc(); iC++) {
    clout << "Cuboid #" << iC << ": " << std::endl;
    get(iC).print();
  }
}

template<typename T>
void CuboidGeometry3D<T>::writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer)
{
  std::ofstream fout;
  if ( singleton::mpi().isMainProcessor() ) {

    // Open File
    fout.open(completeFileName.c_str(), std::ios::app);
    if (!fout) {
      clout << "Error: could not open " << completeFileName << std::endl;
    }

    // --- Preamble --- //
    fout << "<CuboidGeometry dimension=\"3\" " << _cuboidParameters(getMotherCuboid()) << ">\n";

    // TODO: Move Cuboid XML Serialization to Cuboid3D class
    for (int iC = 0; iC < getNc(); ++iC) {
      fout << "<Cuboid " << _cuboidParameters(get(iC)) << " />\n";
    }

    fout << "</CuboidGeometry>\n";

    // Close File
    fout.close();
  }
}


template<typename T>
void CuboidGeometry3D<T>::writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer)
{
  std::string fname = singleton::directories().getLogOutDir() + fileName + ".xml";
  std::ofstream fout;
  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::trunc);
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<XMLContent>\n";
    fout.close();
    fout.clear();
  }

  writeToExistingFile(fname, loadBalancer);

  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::app);
    fout << "</XMLContent>\n";
    fout.close();
  }
}


// TODO: Move this method to Cuboid3D<T> class
/// Helper Function to create cuboid parameters for XML tag
template<typename T>
std::string CuboidGeometry3D<T>::_cuboidParameters(Cuboid3D<T> const& cub)
{
  std::stringstream ss;
  ss.flags(std::ios::scientific);
  ss.precision (std::numeric_limits<double>::digits10 + 1);
  ss << " extent=\"";
  for (int i = 0; i<3; i++) {
    ss << cub.getExtent()[i] << " ";
  }

  ss << "\" origin=\"";
  for (int i = 0; i<3; i++) {
    ss << cub.getOrigin()[i] << " ";
  }

  ss << "\" deltaR=\"" << cub.getDeltaR();
  ss << "\" weight=\"" << cub.getWeightValue();
  ss << "\" refinementLevel=\"" << cub.getRefinementLevel() << "\"";
  return ss.str();
}


}  // namespace olb

#endif
