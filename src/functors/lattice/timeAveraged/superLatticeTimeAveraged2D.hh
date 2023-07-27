/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Mathias J. Krause, Benedict Hasenauer
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

#ifndef SUPER_LATTICE_TIME_AVERAGED_F2_D_HH
#define SUPER_LATTICE_TIME_AVERAGED_F2_D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<limits>
#include "superLatticeTimeAveraged2D.h"


namespace olb {
template <typename T>
SuperLatticeTimeAveragedF2D<T>:: SuperLatticeTimeAveragedF2D( SuperF2D<T,T>& sFunctor)
  : SuperF2D<T,T>(sFunctor.getSuperStructure(),sFunctor.getTargetDim()*2), _ensembles(0), _sFunctor(sFunctor),
    _sData(_sFunctor.getSuperStructure().getCuboidGeometry(),
           _sFunctor.getSuperStructure().getLoadBalancer(),
           _sFunctor.getSuperStructure().getOverlap(),
           _sFunctor.getTargetDim()),
    _sDataP2(_sFunctor.getSuperStructure().getCuboidGeometry(),
             _sFunctor.getSuperStructure().getLoadBalancer(),
             _sFunctor.getSuperStructure().getOverlap(),
             _sFunctor.getTargetDim())
{
  this->getName() = "Time Averaged " + _sFunctor.getName();
};
template <typename T>
bool SuperLatticeTimeAveragedF2D<T>::operator() (T output[], const int input[])
{
  T iCloc = _sData.getLoadBalancer().loc(input[0]);
  for ( int iDim = 0; iDim < _sData.getDataSize(); iDim++) {
    output[iDim] = _sData.getBlock(iCloc).get(input+1,iDim) / _ensembles;
  }
  for (int iDim = _sData.getDataSize(); iDim < _sData.getDataSize()*2; iDim++)
    if (_sDataP2.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())/_ensembles - _sData.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())*_sData.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())/_ensembles/_ensembles<0) {
      output[iDim]=0;
    }
    else {
      output[iDim] = util::sqrt(_sDataP2.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())/_ensembles - _sData.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())*_sData.getBlock(iCloc).get(input+1,(int) iDim-_sDataP2.getDataSize())/_ensembles/_ensembles);
    }
  return true;
};
template <typename T>
int SuperLatticeTimeAveragedF2D<T>::getEnsembles()
{
  return _ensembles;
};
template <typename T>
void SuperLatticeTimeAveragedF2D<T>::addEnsemble()
{
  int i[3];
  for (int iCloc=0; iCloc < _sData.getLoadBalancer().size(); ++iCloc) {
    i[0] = _sData.getLoadBalancer().glob(iCloc);
    _sData.getBlock(iCloc).forSpatialLocations([&](auto iX, auto iY) {
      i[1] = iX;
      i[2] = iY;
      std::vector<BaseType<T>> tmp(_sFunctor.getTargetDim(), 0);
      _sFunctor(tmp.data(), i);
      for (int iDim=0; iDim<_sFunctor.getTargetDim(); iDim++) {
        _sData.getBlock(iCloc).get({iX, iY}, iDim) += (BaseType<T>)(tmp[iDim]) ;
        _sDataP2.getBlock(iCloc).get({iX, iY}, iDim) += (BaseType<T>)(tmp[iDim]) *(BaseType<T>)(tmp[iDim]) ;
      }
    });
  }
  _ensembles++;
};
template <typename T>
int SuperLatticeTimeAveragedF2D<T>::getBlockFSize() const
{
  return 0;
};

template <typename T>
SuperLatticeTimeAveragedCrossCorrelationF2D<T>::SuperLatticeTimeAveragedCrossCorrelationF2D(SuperF2D<T,T>& sFunctorM,SuperF2D<T,T>& sFunctorN)
  : SuperF2D<T,T>(sFunctorM.getSuperStructure(),sFunctorM.getTargetDim()*sFunctorN.getTargetDim()), _ensembles(0), _sFunctorM(sFunctorM), _sFunctorN(sFunctorN), _sDataM(_sFunctorM.getSuperStructure().getCuboidGeometry(),_sFunctorM.getSuperStructure().getLoadBalancer(),_sFunctorM.getSuperStructure().getOverlap(),_sFunctorM.getTargetDim()),_sDataN(_sFunctorN.getSuperStructure().getCuboidGeometry(),_sFunctorN.getSuperStructure().getLoadBalancer(),_sFunctorN.getSuperStructure().getOverlap(),_sFunctorN.getTargetDim()),_sDataMN(_sFunctorM.getSuperStructure().getCuboidGeometry(),_sFunctorM.getSuperStructure().getLoadBalancer(),_sFunctorM.getSuperStructure().getOverlap(),_sFunctorM.getTargetDim()*_sFunctorN.getTargetDim())
{
  this->getName() = "Time Averaged Corss Correlation " + _sFunctorM.getName()+"-"+_sFunctorN.getName();
};

template <typename T>
void SuperLatticeTimeAveragedCrossCorrelationF2D<T>::addEnsemble()
{
  int i[3];
  int iX,iY;
  int iDimMN;


  for (int iCloc=0; iCloc < _sDataMN.getLoadBalancer().size(); ++iCloc) {
    i[0] = _sDataMN.getLoadBalancer().glob(iCloc);
    for (iX=0; iX < _sDataMN.getBlock(iCloc).getNx(); iX++) {
      for (iY=0; iY < _sDataMN.getBlock(iCloc).getNy(); iY++) {
        i[1] = iX;
        i[2] = iY;
        BaseType<T> tmpN[_sFunctorN.getTargetDim()];
        BaseType<T> tmpM[_sFunctorM.getTargetDim()];
        _sFunctorN(tmpN, i);
        _sFunctorM(tmpM, i);
        iDimMN=0;
        for (int iDimM=0; iDimM<_sFunctorM.getTargetDim(); iDimM++) {
          for (int iDimN=0; iDimN<_sFunctorN.getTargetDim(); iDimN++) {
            _sDataMN.getBlock(iCloc).get({iX, iY}, iDimMN) += (BaseType<T>)(tmpM[iDimM])*(BaseType<T>)(tmpN[iDimN]) ;
            iDimMN++;
          }
        }
        for (int iDim=0; iDim<_sFunctorN.getTargetDim(); iDim++) {
          _sDataN.getBlock(iCloc).get({iX, iY}, iDim) += (BaseType<T>)(tmpN[iDim]) ;
        }
        for (int iDim=0; iDim<_sFunctorM.getTargetDim(); iDim++) {
          _sDataM.getBlock(iCloc).get({iX, iY}, iDim) += (BaseType<T>)(tmpM[iDim]) ;
        }
      }
    }
  }

  _ensembles++;
};
template <typename T>
bool SuperLatticeTimeAveragedCrossCorrelationF2D<T>::operator() (T output[], const int input[])
{
  int iDim =0;
  T iCloc = _sDataMN.getLoadBalancer().loc(input[0]);
  for (int iDimM=0; iDimM<_sFunctorM.getTargetDim(); iDimM++) {
    for (int iDimN=0; iDimN<_sFunctorN.getTargetDim(); iDimN++) {
      output[iDim] = _sDataMN.getBlock(iCloc).get(input+1,iDim)-_sDataM.getBlock(iCloc).get(input+1,iDimM) *_sDataN.getBlock(iCloc).get(input+1,iDimN)/_ensembles/_ensembles;
      iDim++;
    }
  }

  return true;

};
template <typename T>
SuperLatticeTimeAveraged2DL2Norm<T>::SuperLatticeTimeAveraged2DL2Norm(SuperF2D<T,T>& sFunctorM,SuperF2D<T,T>& sFunctorN,SuperGeometry<T,2>& sGeometry,int material)
  : SuperF2D<T,T>(sFunctorM.getSuperStructure(),sFunctorM.getTargetDim()), _sFunctorM(sFunctorM), _sFunctorN(sFunctorN), _sGeometry(sGeometry),_material(material)
{
  this->getName() = "SuperLatticeTimeAveraged2DL2Norm";
};

template <typename T>
bool SuperLatticeTimeAveraged2DL2Norm<T>::operator() (T output[], const int input[])
{
  output[0]=0;
  CuboidGeometry2D<T>& geometry = _sFunctorM.getSuperStructure().getCuboidGeometry();

  int inputTmp[3];
  T tmpM[_sFunctorM.getTargetDim()];
  T tmpN[_sFunctorN.getTargetDim()];
  for (int iC = 0; iC < _sFunctorM.getSuperStructure().getLoadBalancer().size(); ++iC) {
    Cuboid2D<T>& cuboid = geometry.get(_sFunctorM.getSuperStructure().getLoadBalancer().glob(iC));

    const int nX     = cuboid.getNx();
    const int nY     = cuboid.getNy();

    inputTmp[0] = _sFunctorM.getSuperStructure().getLoadBalancer().glob(iC);

    for (inputTmp[1] = 0; inputTmp[1] < nX; ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < nY; ++inputTmp[2]) {
        _sFunctorM(tmpM, inputTmp);
        _sFunctorN(tmpN, inputTmp);
        for (int iDim = 0; iDim < _sFunctorM.getTargetDim()/2; ++iDim) {
          output[0]  += (tmpM[iDim]-tmpN[iDim])*(tmpM[iDim]-tmpN[iDim]);
        }

      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0],MPI_SUM);
#endif

  Cuboid2D<T>& cuboid = geometry.get(_sFunctorM.getSuperStructure().getLoadBalancer().glob(0));
  const T   weight = cuboid.getDeltaR();

  output[0]=util::sqrt(output[0])*weight;
  return true;

};
}

#endif
