/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2013, 2016 Mathias J. Krause, Lukas Baron, Benjamin Förster
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


// this file contains the controlled functions that are used to set up the
// fluid flow problem in dependence of a control variable, e.g. the boundary
// velocities, the  pressure and the porosity

#ifndef CONTROLLED_FUNCTIONS_3D_H
#define CONTROLLED_FUNCTIONS_3D_H

#include "utilities/omath.h"
#include <vector>
#include <map>

#include "controller.h"
#include "dualFunctors3D.h"
#include "dualFunctors3D.hh"

#include "functors/functors3D.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {


/// controlled analytical functor
template <typename T, typename DESCRIPTOR>
class ControlledAnalyticalF3D : public AnalyticalF3D<T,T> {
protected:
  Controller<T>* _controller;
  UnitConverter<T,DESCRIPTOR>* _converter;
  SuperLattice<T,DESCRIPTOR>* _lattice;
  SuperGeometry<T,3>* _geometry;


public:
  ControlledAnalyticalF3D(SuperLattice<T,DESCRIPTOR>* lattice, SuperGeometry<T,3>* geometry,
                          Controller<T>* controller, UnitConverter<T,DESCRIPTOR>* converter, int targetDim)
    : AnalyticalF3D<T,T>(targetDim), _controller(controller), _converter(converter),
      _lattice(lattice), _geometry(geometry)
  { }

  virtual std::string name() = 0;

};


/// controlled analytical Velocity function
template <typename T, typename DESCRIPTOR>
class ControlledAnalyticalVelocity3D : public ControlledAnalyticalF3D<T,DESCRIPTOR> {
private:

public:
  ControlledAnalyticalVelocity3D(SuperLattice<T,DESCRIPTOR>* lattice, SuperGeometry<T,3>* geometry, Controller<T>* controller, UnitConverter<T,DESCRIPTOR>* converter)
    : ControlledAnalyticalF3D<T,DESCRIPTOR>(lattice, geometry, controller, converter,3) { }

  bool operator()(T output[], const T input[]);
  std::string name()
  {
    return "startVelocity";
  }

};

/// controlled analytical pressure function
template <typename T, typename DESCRIPTOR>
class ControlledAnalyticalPressure3D : public ControlledAnalyticalF3D<T,DESCRIPTOR> {
private:

public:
  ControlledAnalyticalPressure3D(SuperLattice<T,DESCRIPTOR>* lattice, SuperGeometry<T,3>* geometry, Controller<T>* controller, UnitConverter<T,DESCRIPTOR>* converter)
    : ControlledAnalyticalF3D<T,DESCRIPTOR>(lattice,geometry, controller, converter,1) { }

  bool operator()(T output[], const T input[]);
  std::string name()
  {
    return "startPressure";
  }

};


/// controlled analytical ExternalField
template <typename T, typename DESCRIPTOR>
class ControlledAnalyticalExternalField3D : public ControlledAnalyticalF3D<T,DESCRIPTOR> {
private:

public:
  ControlledAnalyticalExternalField3D(SuperLattice<T,DESCRIPTOR>* lattice, SuperGeometry<T,3>* geometry, Controller<T>* controller, UnitConverter<T,DESCRIPTOR>* converter, int targetDim = 3)
    : ControlledAnalyticalF3D<T,DESCRIPTOR>(lattice,geometry, controller, converter, targetDim) { }

  bool operator()(T output[], const T input[]);
  std::string name()
  {
    return "startExternalField";
  }
};

/// controlled super lattice functor
template <typename T, typename DESCRIPTOR>
class ControlledSuperLatticeF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  Controller<T>& _controller;
  UnitConverter<T,DESCRIPTOR>& _converter;
  SuperGeometry<T,3>& _geometry;

public:
  ControlledSuperLatticeF3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry,
                            Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : SuperLatticeF3D<T,DESCRIPTOR>(lattice,targetDim),
      _controller(controller), _converter(converter), _geometry(geometry)
  { };
};

template <typename T, typename DESCRIPTOR>
class Objective3D : public ControlledSuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperLatticeF3D<T, DESCRIPTOR>* _referenceSolution;
  SuperIndicatorF3D<T>* _objectiveDomain;
public:
  Objective3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry,
              Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter,
              SuperLatticeF3D<T, DESCRIPTOR>* referenceSolution, SuperIndicatorF3D<T>* objectiveDomain)
    : ControlledSuperLatticeF3D<T,DESCRIPTOR>(lattice, geometry, controller, converter, 1),
      _referenceSolution(referenceSolution), _objectiveDomain(objectiveDomain) {
    this->getName() = "objective";
  }

  bool operator()(T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class DObjectiveDf3D : public ControlledSuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperF3D<T,T>* _dObjectiveDf;
  SuperF3D<T,T>* _dObjectiveDf_yzPlane;
  SuperLatticeF3D<T, DESCRIPTOR>* _referenceSolution;
  SuperIndicatorF3D<T>* _objectiveDomain;

  SuperLatticePhysVelocity3D<T,DESCRIPTOR> physLatticeVelocity;
  SuperLatticeDphysVelocityDf3D<T,DESCRIPTOR> dPhysVelocityf;
public:
  DObjectiveDf3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry,
                 Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter,
                 SuperLatticeF3D<T, DESCRIPTOR>* referenceSolution, SuperIndicatorF3D<T>* objectiveDomain);

  ~DObjectiveDf3D() {
    delete _dObjectiveDf;
  }
  bool operator()(T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class DObjectiveDControl3D : public ControlledSuperLatticeF3D<T,DESCRIPTOR> {
protected:

public:
  DObjectiveDControl3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry,
                       Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter,
                       int targetDim)
    : ControlledSuperLatticeF3D<T,DESCRIPTOR>(lattice, geometry, controller, converter, targetDim) {
    this->getName() = "DObjectiveDControl";
  };
  bool operator()(T output[], const int input[]);
};


// ---------------------------
// Define Projections for the control variables
// ---------------------------

/// Projection Base Class
/**
 * A Projection maps an alpha \in (-\infty, \infty) on a control \in [0,1].
 */
template <typename T, typename DESCRIPTOR>
class Projection3D : public ControlledSuperLatticeF3D<T,DESCRIPTOR> {

protected:
  int _material;
  int _nVoxel;
  std::vector<int> _nVoxelCuboid;
  int _nCuboids;
  int _nDim;
  T _scale;

public:
  Projection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
               Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim, T scale=T(1))
    : ControlledSuperLatticeF3D<T,DESCRIPTOR>(lattice, geometry, controller, converter, targetDim),
      _material(material), _scale(scale)
  {
    this->getName() = "alpha";
    _nVoxel = geometry.getStatistics().getNvoxel(material);
    _nCuboids = geometry.getCuboidGeometry().getNc();
    _nVoxelCuboid.resize(_nCuboids);
    _nVoxelCuboid[0] = int();
    for (int iC=1; iC<_nCuboids; iC++) {
      _nVoxelCuboid[iC] = _nVoxelCuboid[iC-1] + geometry.getCuboidGeometry().get(iC-1).getLatticeVolume();
    }
    _nDim = targetDim;
  }

  /// Copy Constructor
  Projection3D(Projection3D<T,DESCRIPTOR>& rhs)
    : Projection3D (rhs._sLattice, rhs._geometry, rhs._material, rhs._controller, rhs._converter, rhs._nDim, rhs._scale)
  { }

  virtual T getInverseFunction(T porosity) {
    std::cout << "WARNING: no inverse for projection found! Exiting..." << std::endl;
    exit(-1);
  }

  int getIndex(const int input[], int iDim)
  {
    int nX = this->_geometry.getCuboidGeometry().get(input[0]).getNx();
    int nY = this->_geometry.getCuboidGeometry().get(input[0]).getNy();
    /// [NDIM*NX*NY*NZ]*c + [NDIM*NX*NY]*z + [NDIM*NX]*y + [NDIM]*x + dim
    return _nDim*_nVoxelCuboid[input[0]] + _nDim*nX*nY*input[3] + _nDim*nX*input[2] + _nDim*input[1] + iDim;
  }

  T& getControl(const int input[], int iDim) {
    return this->_controller.getControl(getIndex(input, iDim) );
  }

  T gridTerm() {
    return util::pow( this->_converter.getPhysDeltaX(), 2) * this->_converter.getLatticeViscosity() * this->_converter.getLatticeRelaxationTime();
  }

  bool operator()(T output[], const int input[]) {
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      output[iDim] = getControl(input, iDim) * _scale;
    }
    return true;
  }
};




// ------------------------------------------ //
//       grid-independent projection         //
// ------------------------------------------ //


/// Projection based on the logistic sigmoid function
template <typename T, typename DESCRIPTOR>
class SigmoidFunction : public Projection3D<T,DESCRIPTOR> {
public:
  SigmoidFunction(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                  Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl; // alpha
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = util::exp(ctrl) / ( util::exp(ctrl) + 1 );
    }
    return true;
  };

  /// Inverse sigmoid function to get control startValue from porosity
  T getInverseFunction(T porosity)
  {
    return -util::log(1./porosity - 1.);
  }
};


/// Derivative of sigmoid function with respect to alpha
template <typename T, typename DESCRIPTOR>
class DSigmoidFunctionDAlpha : public Projection3D<T,DESCRIPTOR> {
public:
  DSigmoidFunctionDAlpha(Projection3D<T,DESCRIPTOR>& proj)
    : Projection3D<T,DESCRIPTOR>(proj)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = util::exp(ctrl) / util::pow( util::exp(ctrl) + 1, 2. );
    }
    return true;
  };
};


/// Projection based on the rectifier function
template <typename T, typename DESCRIPTOR>
class RectifierFunction : public Projection3D<T,DESCRIPTOR> {
public:
  RectifierFunction(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                    Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl; // alpha
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      //ctrl = util::max( this->gridTerm(), ctrl );
      //output[iDim] = 1. - this->gridTerm()/ctrl;
      ctrl = util::max( 1., ctrl );
      // output will be the porosity (startExternalField)
      // using the grid term this corresponds to d=1-G_h/K
      output[iDim] = 1. - 1./ctrl;
    }
    return true;
  };

//  /// Inverse sigmoid function to get control startValue from porosity
//  T getInverseFunction(T porosity) { }
};


/// Derivative of rectifier function with respect to alpha
template <typename T, typename DESCRIPTOR>
class DRectifierFunctionDAlpha : public Projection3D<T,DESCRIPTOR> {
public:
  DRectifierFunctionDAlpha(Projection3D<T,DESCRIPTOR>& proj)
    : Projection3D<T,DESCRIPTOR>(proj)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      //if (ctrl > this->gridTerm()) {
      if (ctrl > 1.) {
        output[iDim] = 1.;
      }
      else {
        output[iDim] = 0.;
      }
    }
    return true;
  };
};



/// Projection based on smooth approximation of rectifier function (Softplus)
template <typename T, typename DESCRIPTOR>
class SoftplusFunction : public Projection3D<T,DESCRIPTOR> {
public:
  SoftplusFunction(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                   Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl; // alpha
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      //output[iDim] = util::log( 1 + util::exp( ctrl - 50 ) ) + this->gridTerm();
      ctrl = util::log( 1 + util::exp( ctrl )) + 1.;
      // output will be the porosity (startExternalField)
      // using the grid term this corresponds to d=1-G_h/K
      output[iDim] = 1. - 1./ctrl;
    }
    return true;
  };

//  /// Inverse sigmoid function to get control startValue from porosity
//  T getInverseFunction(T porosity) { }
};


/// Derivative of softplus function with respect to alpha
template <typename T, typename DESCRIPTOR>
class DSoftplusFunctionDAlpha : public Projection3D<T,DESCRIPTOR> {
public:
  DSoftplusFunctionDAlpha(Projection3D<T,DESCRIPTOR>& proj)
    : Projection3D<T,DESCRIPTOR>(proj)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      //output[iDim] = util::exp(ctrl) / ( util::exp(ctrl) + util::exp(50.) );
      output[iDim] = util::exp(ctrl) / ( util::exp(ctrl) + 1. );
    }
    return true;
  };
};


/// Lukas Baron Projection
/**
 * Projection Mapping: d(a) = (0.5 + 0.5 * util::sin ( ( a - 0.5 ) * PI ) ** 2
 */
template <typename T, typename DESCRIPTOR>
class BaronProjection3D : public Projection3D<T,DESCRIPTOR> {
public:
  BaronProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                    Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      output[iDim] = util::pow( .5 + .5 * util::sin(
                                  ( this->_controller.getControl( this->getIndex(input, iDim) ) - .5 ) * 3.14 ), 2.);
    }
    return true;
  };
};


/// Derivative of BaronProjection3D with respect to alpha
template <typename T, typename DESCRIPTOR>
class DBaronProjectionDAlpha3D : public Projection3D<T,DESCRIPTOR> {

public:
  DBaronProjectionDAlpha3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                           Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };
  bool operator()(T output[], const int input[])
  {
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      output[iDim] = ( .5 + .5 * util::sin(
                         ( this->_controller.getControl( this->getIndex(input, iDim) ) - .5) * 3.14 ) )
                     * util::cos( ( this->_controller.getControl( this->getIndex(input, iDim) ) - .5 ) * 3.14 )
                     * 3.14;
    }
    return true;
  };
};


/// Mathias Krause Projection
/**
 * Projection Mapping: d(a) = util::exp(a**2)
 */
template <typename T, typename DESCRIPTOR>
class KrauseProjection3D : public Projection3D<T,DESCRIPTOR> {
public:
  KrauseProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                     Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = 1 - ( 1. / util::pow( 2.7182, util::pow( ctrl, 2. ) ) );
    }
    return true;
  };
};


/// Derivative of KrauseProjection3D with respect to alpha
template <typename T, typename DESCRIPTOR>
class DKrauseProjectionDAlpha3D : public Projection3D<T,DESCRIPTOR> {

public:
  DKrauseProjectionDAlpha3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                            Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {

    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = ( 2 * ctrl / util::pow( 2.7182, util::pow( ctrl, 2. ) ) );
    }
    return true;
  };
};


/********* DEPRECATED (see new GIProjections) ************/
//template <typename T, typename DESCRIPTOR>
//class FoersterProjection3D : public Projection3D<T,DESCRIPTOR> {
//public:
//  FoersterProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
//                 Controller<T>& controller, LBconverter<T>& converter, int targetDim)
//      : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
//  { };
//
//  bool operator()(T output[], const int input[]) {
//
//    T ctrl; // alpha
//    //T K; // permeability
//    for (int iDim = 0; iDim < this->_nDim; iDim++) {
//      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
//      //K = util::pow( 2.71828182, ctrl) + tau * nuLattice * util::pow( physLength, T( dim-1 ) ) ;
//      output[iDim] = util::exp(ctrl) / ( util::exp(ctrl) + this->gridTerm() );
//    }
//    return true;
//  };
//};
//
//template <typename T, typename DESCRIPTOR>
//class DFoersterProjectionDControl3D : public Projection3D<T,DESCRIPTOR> {
//
//public:
//  DFoersterProjectionDControl3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
//                          Controller<T>& controller, LBconverter<T>& converter, int targetDim)
//      : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
//  { };
//
//  bool operator()(T output[], const int input[]) {
//
//    T ctrl;
//
//    for (int iDim = 0; iDim < this->_nDim; iDim++) {
//      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
//      output[iDim] = util::exp(ctrl) * this->gridTerm() / util::pow( util::exp(ctrl) + this->gridTerm(), 2. );
//    }
//    return true;
//  };
//};
//
//
//
//template <typename T, typename DESCRIPTOR>
//class Foerster2Projection3D : public Projection3D<T,DESCRIPTOR> {
//public:
//  Foerster2Projection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
//                       Controller<T>& controller, LBconverter<T>& converter, int targetDim)
//      : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
//  { };
//
//  bool operator()(T output[], const int input[]) {
//
//    T ctrl; // alpha
//    //T K; // permeability
//    for (int iDim = 0; iDim < this->_nDim; iDim++) {
//      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
//      //K = util::pow( 2.71828182, ctrl) + tau * nuLattice * util::pow( physLength, T( dim-1 ) ) ;
//      output[iDim] = (util::exp(util::pow(ctrl, 2.0)) - 1.0) / ( util::exp(util::pow(ctrl, 2.0)) - 1.0 + this->gridTerm() );
//    }
//    return true;
//  };
//};
//
//template <typename T, typename DESCRIPTOR>
//class DFoerster2ProjectionDControl3D : public Projection3D<T,DESCRIPTOR> {
//
//public:
//  DFoerster2ProjectionDControl3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
//                                Controller<T>& controller, LBconverter<T>& converter, int targetDim)
//      : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
//  { };
//
//  bool operator()(T output[], const int input[]) {
//
//    T ctrl;
//
//    for (int iDim = 0; iDim < this->_nDim; iDim++) {
//      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
//      output[iDim] = 2* ctrl * util::exp(util::pow(ctrl, 2.0)) * this->gridTerm() / util::pow( util::exp(util::pow(ctrl, 2.0)) - 1.0 +
//                                                                                                this->gridTerm(), 2. );
//    }
//    return true;
//  };
//};



/// Simon Stasius Projection (DEPRECATED! See new GIProjections)
template <typename T, typename DESCRIPTOR>
class StasiusProjection3D : public Projection3D<T,DESCRIPTOR> {
public:
  StasiusProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                      Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T nuLattice = this->_converter.getLatticeNu();
    T tau =  this->_converter.getTau();
    T physLength = this->_converter.physLength();
    int dim = this->_nDim; //this->_converter.getDim();

    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = util::pow( ctrl, 2. ) / ( util::pow( ctrl, 2. ) +  tau * nuLattice * util::pow( physLength, T( dim-1 ) ) );
    }
    return true;
  };
};


/// Simon Stasius Projection (DEPRECATED! See new GIProjections)
template <typename T, typename DESCRIPTOR>
class DStasiusProjectionDAlpha3D : public Projection3D<T,DESCRIPTOR> {

public:
  DStasiusProjectionDAlpha3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                             Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  bool operator()(T output[], const int input[])
  {
    T nuLattice = this->_converter.getLatticeNu();
    T tau =  this->_converter.getTau();
    T physLength = this->_converter.physLength();
    int dim = 3; //this->_converter.getDim();

    T const_term = tau * nuLattice * util::pow( physLength, T( dim-1 ) );

    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = 2. * ctrl * const_term / util::pow( ( util::pow( ctrl, 2. ) + const_term ), 2.);
    }
    return true;
  };
};


/// Gridterm-dependent Projection Base Class
/**
 * Ba = 1 - GridTerm() / K(a)
 * with K(a) = subprojection(a) + GridTerm()
 */
template <typename T, typename DESCRIPTOR>
class GiProjection3D : public Projection3D<T,DESCRIPTOR> {
public:
  GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                 Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : Projection3D<T,DESCRIPTOR>(lattice, geometry, material, controller, converter, targetDim, 1)
  { };

  virtual T subprojection(T ctrl) = 0;
  virtual T dSubprojection(T ctrl) = 0;

  virtual bool operator()(T output[], const int input[])
  {
    T ctrl; // alpha
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = subprojection(ctrl) / ( subprojection(ctrl) + this->gridTerm() );
    }
    return true;
  };
};


/// Derivative of GiProjection3D with respect to alpha
template <typename T, typename DESCRIPTOR>
class DGiProjectionDAlpha3D : public Projection3D<T,DESCRIPTOR> {
protected:
  GiProjection3D<T,DESCRIPTOR>& _giProjection;
public:
  DGiProjectionDAlpha3D(GiProjection3D<T,DESCRIPTOR>& proj)
    : Projection3D<T,DESCRIPTOR>(proj),
      _giProjection(proj)
  { };

  virtual bool operator()(T output[], const int input[])
  {
    T ctrl;
    for (int iDim = 0; iDim < this->_nDim; iDim++) {
      ctrl = this->_controller.getControl( this->getIndex(input, iDim) );
      output[iDim] = _giProjection.dSubprojection(ctrl) * this->gridTerm() /
                     util::pow( _giProjection.subprojection(ctrl) + this->gridTerm(), 2. );
    }
    return true;
  };
};


template <typename T, typename DESCRIPTOR>
class FoersterGiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  FoersterGiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                         Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::exp(ctrl);
  }

  T dSubprojection(T ctrl)
  {
    return util::exp(ctrl);
  }
};


/// FoersterProjection for arbitrary n
/**
 * subproj(a) = util::exp(a^(2n)) - 1
 */
template <typename T, typename DESCRIPTOR>
class FoersterNGiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
protected:
  int _n; /// n in util::exp(a^(2n)) - 1
public:
  FoersterNGiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                          Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int n, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim),
      _n(n)
  {
    if (n < 1) {
      std::cout << "Error! n must be at least 1. Exiting..." << std::endl;
      exit(-1);
    }
  };
  T subprojection(T ctrl)
  {
    return util::exp( util::pow(ctrl, 2.0 * _n) ) - 1.0;
  }

  T dSubprojection(T ctrl)
  {
    return 2.0*_n * util::pow(ctrl, (2.0*_n - 1.0)) * util::exp( util::pow(ctrl, 2.0 * _n) );
  }
};


/// StasiusProjection for arbitrary n
/**
 * subproj(a) = a^(2n)
 */
template <typename T, typename DESCRIPTOR>
class StasiusNGiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
protected:
  int _n; /// n in util::exp(a^(2n))
public:
  StasiusNGiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                         Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int n, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim),
      _n(n)
  {
    if (n < 1) {
      std::cout << "Error! n must be at least 1. Exiting..." << std::endl;
      exit(-1);
    }
  };
  T subprojection(T ctrl)
  {
    return util::pow(ctrl, 2.0 * _n);
  }

  T dSubprojection(T ctrl)
  {
    return 2.0*_n * util::pow(ctrl, (2.0*_n - 1.0));
  }
};


template <typename T, typename DESCRIPTOR>
class Foerster2GiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  Foerster2GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                          Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::exp( util::pow(ctrl, 2.0) ) - 1.0;
  }

  T dSubprojection(T ctrl)
  {
    return 2.0 * ctrl * util::exp( util::pow(ctrl, 2.0) );
  }
};


template <typename T, typename DESCRIPTOR>
class Foerster3GiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  Foerster3GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                          Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::exp( util::pow(ctrl, 4.0) ) - 1.0;
  }

  T dSubprojection(T ctrl)
  {
    return 4.0 * util::pow(ctrl, 3.0) * util::exp( util::pow(ctrl, 4.0) );
  }
};


template <typename T, typename DESCRIPTOR>
class Foerster4GiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  Foerster4GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                          Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::exp( util::pow(ctrl, 10.0) ) - 1.0;
  }

  T dSubprojection(T ctrl)
  {
    return 10.0 * util::pow(ctrl, 9.0) * util::exp( util::pow(ctrl, 10.0) );
  }
};


template <typename T, typename DESCRIPTOR>
class StasiusGiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  StasiusGiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                        Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::pow(ctrl, 2.0);
  }

  T dSubprojection(T ctrl)
  {
    return 2.0 * ctrl;
  }
};



template <typename T, typename DESCRIPTOR>
class Stasius2GiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  Stasius2GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                         Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::pow(ctrl, 4.0);
  }

  T dSubprojection(T ctrl)
  {
    return 4.0 * util::pow(ctrl, 3.0);
  }
};


template <typename T, typename DESCRIPTOR>
class Stasius3GiProjection3D : public GiProjection3D<T,DESCRIPTOR> {
public:
  Stasius3GiProjection3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry, int material,
                         Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
    : GiProjection3D<T,DESCRIPTOR> (lattice, geometry,  material, controller, converter, targetDim)
  { };
  T subprojection(T ctrl)
  {
    return util::pow(ctrl, 10.0);
  }

  T dSubprojection(T ctrl)
  {
    return 10.0 * util::pow(ctrl, 9.0);
  }
};




/*
////////////////////////////////////////////////////////////////////////////////
/// functor to compute the objective functional
template<typename T, typename DESCRIPTOR>
class ComputeVolumeFlowPenalty : public LatticeF3D<T,DESCRIPTOR> {
private:
  SuperLattice<T, DESCRIPTOR>* sLattice;
  SuperGeometry<T,3>* bg;
  Dynamics<T, DESCRIPTOR>* bulkDynamics;
  LBconverter<T>* converter;
  T u0;
  T alpha;

public:
  // Constructor
  ComputeVolumeFlowPenalty(SuperLattice<T, DESCRIPTOR>* _sLattice,
                      SuperGeometry<T,3>* _bg,
                      Dynamics<T, DESCRIPTOR>* _bulkDynamics,
                      LBconverter<T>* _converter,
                      T _u0,
                      T _alpha
                      )
    : LatticeF3D<T,DESCRIPTOR> (NULL,NULL),
      sLattice(_sLattice),
      bg(_bg),
      bulkDynamics(_bulkDynamics),
      converter(_converter),
      u0(_u0), alpha(_alpha)
  { }

  std::vector<T> operator()(std::vector<T> input);

  // obligatory access operator


};*/

template <typename T, typename DESCRIPTOR>
class RelativeObjective3D : public Objective3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>& _solution;
  SuperIndicatorF3D<T>& _objectiveDomain;

public:
  RelativeObjective3D(SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& geometry,
                      Controller<T>& controller, UnitConverter<T,DESCRIPTOR>& converter, AnalyticalF3D<T,T>& solution, SuperIndicatorF3D<T>& objectiveDomain)
    : Objective3D<T,DESCRIPTOR>(lattice, geometry, controller, converter), _solution(solution), _objectiveDomain(objectiveDomain)
  {
    this->getName() = "objective";
  };
  bool operator()(T output[], const int input[]);
};

} // namespace opti

} // namespace olb
#endif // CONTROLLED_FUNCTIONS_H

