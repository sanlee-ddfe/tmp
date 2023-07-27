/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster, Julius Jeßberger
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


#ifndef OPTI_CASE_DUAL_HH
#define OPTI_CASE_DUAL_HH

#include "optiCaseDual.h"

#include <regex>

namespace olb {

namespace opti {

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::readFromXML(XMLreader const& xml)
{
  std::string type ("");
  xml.readOrWarn<std::string>("Optimization", "ControlType", "", type);
  _controlType = ((type == "Force") || (type == "force")) ? ForceControl : PorosityControl;
  xml.readOrWarn<std::string>("Optimization", "Projection", "", _projectionName);
  xml.readOrWarn<int>("Optimization", "ControlMaterial", "", _controlMaterial);
  xml.readOrWarn<S>("Optimization", "RegAlpha", "", _regAlpha);
  type = "";
  xml.readOrWarn<std::string>("Optimization", "StartValueType", "", type);
  if (type == "Porosity") {
    _startValueType = Porosity;
  } else if (type == "Permeability") {
    _startValueType = Permeability;
  } else {
    _startValueType = Control;
  }
   xml.readOrWarn<bool>("Optimization", "ReferenceSolution", "", _computeReference);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initialize(XMLreader const& xml)
{
  _converter = createUnitConverter<S,descriptor>(xml);

  _referenceSolver = createLbSolver <SOLVER<S,SolverMode::Reference>> (xml);
  if (_computeReference) {
    _referenceSolver->solve();
  } else {
    // only build up geometry and lattice
    _referenceSolver->buildAndReturn();
  }

  _referenceGeometry = _referenceSolver->parameters(names::Results()).geometry;
  _referenceLattice = _referenceSolver->parameters(names::Results()).lattice;

  if (_controlType == ForceControl) {
    _fieldDim = dim;
  } else {
    _fieldDim = 1;
  }
  _dimCtrl = _referenceGeometry->getStatistics().getNvoxel() * _fieldDim;
  // for mysterious reasons, material 0 is excluded in the function above, so we add it by hand
  _dimCtrl += _referenceGeometry->getStatistics().getNvoxel(0) * _fieldDim;
  _controller = new Controller<S>(_dimCtrl);

  _primalSolver = createLbSolver <SOLVER<S,SolverMode::Primal>> (xml);
  _dualSolver = createLbSolver <SOLVER<S,SolverMode::Dual>> (xml);

  if (_computeReference) {
    _primalSolver->parameters(names::Opti()).referenceSolution
     = std::make_shared<SuperLatticePhysVelocity3D<S,descriptor>>(*_referenceLattice, *_converter);
  }
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::initializeFields()
{
  // set _projection, _dProjectionDcontrol and controlledField of the solvers
  // _projection, _dProjectionDcontrol need the geometry + lattice data, but
  // they are independent of the populations.
  if (_controlType == ForceControl) {
    if (_projectionName == "Sigmoid") {
      _projection = new SigmoidFunction<S, descriptor>(*_referenceLattice, *_referenceGeometry, _controlMaterial,
          *_controller, *_converter, _fieldDim);
      _dProjectionDcontrol = new DSigmoidFunctionDAlpha<S, descriptor>(*_projection);
    }
    else {
      S scale = 1. / (_converter->getConversionFactorForce() / _converter->getConversionFactorMass());
      _projection = new Projection3D<S, descriptor>(*_referenceLattice, *_referenceGeometry, _controlMaterial,
          *_controller, *_converter, _fieldDim, scale);
      AnalyticalConst3D<S,S>* dProjectionDcontrol = new AnalyticalConst3D<S,S> (scale, scale, scale);
      _dProjectionDcontrol = new SuperLatticeFfromAnalyticalF3D<S, descriptor>(dProjectionDcontrol, *_referenceLattice);
    }
  }

  else if (_controlType == PorosityControl) {
    std::smatch match;

    // --- grid-independent projections --- //
    if (_projectionName == "Sigmoid") {
      Projection3D<S, descriptor>* proj_ptr = new SigmoidFunction<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DSigmoidFunctionDAlpha<S, descriptor>(*proj_ptr);
    }
    else if (_projectionName == "Rectifier") {
      Projection3D<S, descriptor>* proj_ptr = new RectifierFunction<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DRectifierFunctionDAlpha<S, descriptor>(*proj_ptr);
    }
    else if (_projectionName == "Softplus") {
      Projection3D<S, descriptor>* proj_ptr = new SoftplusFunction<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DSoftplusFunctionDAlpha<S, descriptor>(*proj_ptr);
    }
    else if (_projectionName == "Baron") {
      _projection = new BaronProjection3D<S, descriptor>(*_referenceLattice, *_referenceGeometry,
          _controlMaterial, *_controller, *_converter, _fieldDim);
      _dProjectionDcontrol = new DBaronProjectionDAlpha3D<S, descriptor>(*_referenceLattice,
          *_referenceGeometry, _controlMaterial,
          *_controller, *_converter, _fieldDim);
    }
    else if (_projectionName == "Krause") {
      _projection = new KrauseProjection3D<S, descriptor>(*_referenceLattice, *_referenceGeometry,
          _controlMaterial, *_controller, *_converter, _fieldDim);
      _dProjectionDcontrol = new DKrauseProjectionDAlpha3D<S, descriptor>(*_referenceLattice,
          *_referenceGeometry, _controlMaterial,
          *_controller, *_converter, _fieldDim);
    }

    // --- gridterm-dependent projections --- //
    else if (_projectionName == "Foerster") {
      GiProjection3D<S, descriptor>* proj_ptr = new FoersterGiProjection3D<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DGiProjectionDAlpha3D<S, descriptor>(*proj_ptr);
    }
    else if (std::regex_match(_projectionName, match, std::regex ("(Foerster)([0-9])"))) {
      int n = std::stoi(match[2]);
      clout << "Creating Projection Foerster" << match[2] << ": util::exp(a^" << 2.0*n << ")" << std::endl;
      GiProjection3D<S, descriptor>* proj_ptr = new FoersterNGiProjection3D<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          n, _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DGiProjectionDAlpha3D<S, descriptor>(*proj_ptr);
    }
    else if (std::regex_match(_projectionName, match, std::regex ("(Stasius)([0-9])"))) {
      int n = std::stoi(match[2]);
      clout << "Creating Projection Stasius" << match[2] << ": a^" << 2.0*n << "" << std::endl;
      GiProjection3D<S, descriptor>* proj_ptr = new StasiusNGiProjection3D<S, descriptor>(*_referenceLattice,
          *_referenceGeometry,
          _controlMaterial, *_controller,
          *_converter,
          n, _fieldDim);
      _projection = proj_ptr;
      _dProjectionDcontrol = new DGiProjectionDAlpha3D<S, descriptor>(*proj_ptr);
    }
  }

  std::shared_ptr<AnalyticalF<dim,S,S>> _analyticalAlpha
   = std::make_shared<SequentialAnalyticalFfromSuperF3D<S>> (*_projection);

  _primalSolver->parameters(names::Opti()).controlledField = _analyticalAlpha;
  _dualSolver->parameters(names::Opti()).controlledField = _analyticalAlpha;
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
S OptiCaseDual<S,SOLVER,C>::evaluateObjective(
  const C& control, unsigned optiStep)
{
  //project control to enforce volume constraint
#ifdef mathiasProjection

  const int nx = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNx();
  const int ny = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNy();
  const int nz = _dualSolver->getGeometry()->getCuboidGeometry().getMotherCuboid().getNz();

  S pi = 4.0 * util::atan(1.0);

  S totalD = T();
  S wantedTotalD = T();

  S upperBound, lowerBound, volumeRatio;
  xml.readOrWarn<T>("Optimization", "UpperBound", "", upperBound);
  xml.readOrWarn<T>("Optimization", "LowerBound", "", lowerBound);
  xml.readOrWarn<T>("Optimization", "VolumeRatio", "", volumeRatio);

  for (int iX=1; iX<nx-1; iX++) {
    for (int iY=1; iY<ny-1; iY++) {
      for (int iZ=1; iZ<nz-1; iZ++) {
        // TODO: this seems to be fixed to material number = 6
        if (_dualSolver->getGeometry()->getMaterial(iX,iY,iZ)==6) {
          S ctrl = control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)];
          totalD += 1./(upperBound-lowerBound)*(ctrl-1.);
          wantedTotalD += 1.0;
        }
      }
    }
  }

  for (int iX=1; iX<nx-1; iX++) {
    for (int iY=1; iY<ny-1; iY++) {
      for (int iZ=1; iZ<nz-1; iZ++) {
        // TODO: this seems to be fixed to material number = 6
        if (_dualSolver->getGeometry()->getMaterial(iX,iY,iZ)==6) {
          control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)]
           = lowerBound
             + (control[(iX-1)*(ny-1)*(nz-1)+(iY-1)*(nz-1)+(iZ-1)] - lowerBound)
               * wantedTotalD * (1.-volumeRatio) / totalD;
        }
      }
    }
  }
#endif

  _controller->setControl(control, _dimCtrl);
  _primalSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;
  _primalSolver->solve();

  return _primalSolver->parameters(names::Results()).objective;
}


template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::computeDerivatives(
  const C& control, C& derivatives, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _dualSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;

  const auto& primalResults = _primalSolver->parameters(names::Results());
  auto& dualParams = _dualSolver->parameters(names::Opti());

  dualParams.fpop = primalResults.fpop;
  dualParams.dObjectiveDf
   = std::make_shared <SequentialAnalyticalFfromSuperF3D<S>> (*primalResults.djdf);
  dualParams.dObjectiveDcontrol
   = std::make_shared <SequentialAnalyticalFfromSuperF3D<S>> (*primalResults.djdalpha);

  _dualSolver->solve();

  derivativesFromDualSolution(derivatives);
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
void OptiCaseDual<S,SOLVER,C>::derivativesFromDualSolution(
  C& derivatives)
{
  const auto& primalGeometry = _primalSolver->parameters(names::Results()).geometry;

  if (_controlType == ForceControl) {
    // get lattice coordinates
    const int nC = primalGeometry->getCuboidGeometry().getNc();
    int latticeR[4];
    for (int iC=0; iC<nC; iC++) {
      latticeR[0] = iC;
      int nX = primalGeometry->getCuboidGeometry().get(iC).getNx();
      int nY = primalGeometry->getCuboidGeometry().get(iC).getNy();
      int nZ = primalGeometry->getCuboidGeometry().get(iC).getNz();
      for (int iX=0; iX<nX; iX++) {
        latticeR[1] = iX;
        for (int iY=0; iY<nY; iY++) {
          latticeR[2] = iY;
          for (int iZ=0; iZ<nZ; iZ++) {
            latticeR[3] = iZ;
            const LatticeR<dim+1> latticeRv_new(latticeR);

            if (getMaterialGlobally(latticeRv_new) == _controlMaterial) {
              C derivativesHelp(_fieldDim, 0);
              if (primalGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
                S rho = _primalSolver->parameters(names::Results()).lattice->get(latticeRv_new).computeRho();
                S dObjectiveDcontrol[_fieldDim];
                (*_primalSolver->parameters(names::Results()).djdalpha)(dObjectiveDcontrol,latticeR);
                S dProjectionDcontrol[_fieldDim];
                (*_dProjectionDcontrol)(dProjectionDcontrol,latticeR);
                for (int iDim=0; iDim<_fieldDim; iDim++) {
                  for (int jPop=0; jPop < descriptor::q; ++jPop) {
                    S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeRv_new)[jPop];
                    derivativesHelp[iDim] -= rho * descriptors::t<S,descriptor>(jPop)
                       * descriptors::invCs2<S,descriptor>() * descriptors::c<descriptor>(jPop,iDim) * phi_j;
                  }
                  // dObjectiveDcontrol is zero if objective does not depend on control
                  // --> dObjectiveDcontrol is 0 and therefore dProjectionDcontrol may be
                  // irrelevant for force and porosity optimisation
                  derivativesHelp[iDim] += _regAlpha * _projection->getControl(latticeR, iDim)
                       + dObjectiveDcontrol[iDim] * dProjectionDcontrol[iDim];
                }
              }
#ifdef PARALLEL_MODE_MPI
              singleton::mpi().bCast(&derivativesHelp[0], _fieldDim, primalGeometry->getLoadBalancer().rank(iC));
#endif
              for (int iDim=0; iDim<_fieldDim; iDim++) {
                derivatives[_projection->getIndex(latticeR,iDim)] = derivativesHelp[iDim];
              }
            }
          }
        }
      }
    }
  }

  if ( _controlType == PorosityControl ) {
    const S omega = _converter->getLatticeRelaxationFrequency();
    const int nC = primalGeometry->getCuboidGeometry().getNc();

    for (int iC = 0; iC < nC; iC++) {

      const Vector<int,3> extend = primalGeometry->getCuboidGeometry().get(iC).getExtent();

      for (int iX = 0; iX < extend[0]; iX++) {
        for (int iY = 0; iY < extend[1]; iY++) {
          for (int iZ = 0; iZ < extend[2]; iZ++) {
            const std::vector<int> latticeR {iC, iX, iY, iZ};
            const LatticeR<dim+1> latticeR_new(iC, iX, iY, iZ);

            if (getMaterialGlobally(latticeR_new) == _controlMaterial) {
              S derivativeHelp (0);

              if (primalGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
                S rho_f = _primalSolver->parameters(names::Results()).lattice->get(latticeR_new).computeRho();
                S u_f[3];
                _primalSolver->parameters(names::Results()).lattice->get(latticeR_new).computeU(u_f);
                S d = _dualSolver->parameters(names::Results()).lattice->get(latticeR_new).template getField<descriptors::POROSITY>();
                //S dObjectiveDcontrol;
                //(*_primalSolver->parameters(names::Results()).dObjectiveDcontrol)(&dObjectiveDcontrol, &latticeR[0]);
                S dProjectionDcontrol;
                (*_dProjectionDcontrol)(&dProjectionDcontrol, &latticeR[0]);

                for (int jPop = 0; jPop < descriptor::q; ++jPop) {
                  S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR_new)[jPop];
                  S feq_j = equilibrium<descriptor>::secondOrder(jPop, rho_f, u_f) + descriptors::t<S,descriptor>(jPop);
                  for (int iDim = 0; iDim < descriptor::d; iDim++) {
                    derivativeHelp +=  phi_j*feq_j*( descriptors::c<descriptor>(jPop,iDim) - d*u_f[iDim] )*u_f[iDim]*dProjectionDcontrol;
                  }
                }
                derivativeHelp *= -omega*descriptors::invCs2<S,descriptor>();
              }

#ifdef PARALLEL_MODE_MPI
              singleton::mpi().bCast(&derivativeHelp, 1, primalGeometry->getLoadBalancer().rank(iC));
#endif
              derivatives[_projection->getIndex(&latticeR[0],0)] = derivativeHelp;
            }
          }
        }
      }
    }
  }
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
C OptiCaseDual<S,SOLVER,C>::getReferenceControl() const
{
  assert((this->_computeReference) && "Reference solution has not been computed\n");
  C result(_dimCtrl, 0);

  if (_controlType == ForceControl) {
    int latticeR[4];
    for (int iC=0; iC<_referenceGeometry->getCuboidGeometry().getNc(); iC++) {
      latticeR[0] = iC;
      int nX = _referenceGeometry->getCuboidGeometry().get(iC).getNx();
      int nY = _referenceGeometry->getCuboidGeometry().get(iC).getNy();
      int nZ = _referenceGeometry->getCuboidGeometry().get(iC).getNz();
      for (int iX=0; iX<nX; iX++) {
        latticeR[1] = iX;
        for (int iY=0; iY<nY; iY++) {
          latticeR[2] = iY;
          for (int iZ=0; iZ<nZ; iZ++) {
            latticeR[3] = iZ;
            const LatticeR<dim+1> latticeRv_new(latticeR);

            if (getMaterialGlobally(latticeRv_new) == _controlMaterial) {
              C force_help(_fieldDim, 0);
              if (_referenceGeometry->getLoadBalancer().rank(iC)  == singleton::mpi().getRank()) {
                for (int iDim=0; iDim<_fieldDim; iDim++) {
                  force_help[iDim]
                   = _referenceLattice->get(latticeRv_new).template getFieldComponent<descriptors::FORCE>(iDim)
                    * _converter->getConversionFactorForce() / _converter->getConversionFactorMass();
                }
              }
#ifdef PARALLEL_MODE_MPI
              singleton::mpi().bCast(&force_help[0], _fieldDim, _referenceGeometry->getLoadBalancer().rank(iC));
#endif
              for (int iDim=0; iDim<_fieldDim; iDim++) {
                result[_projection->getIndex(latticeR,iDim)] = force_help[iDim];
              }
            }
          }
        }
      }
    }
  } else {
    clout << "getReferenceControl() is only implemented for forced optimization" << std::endl;
    exit(1);
  }
  return result;
}

template<typename S, template<typename,SolverMode> typename SOLVER, typename C>
int OptiCaseDual<S,SOLVER,C>::getMaterialGlobally (
  LatticeR<dim+1> latticeR) const
{
  int material = 0;
  if ( _referenceGeometry->getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank() ) {
    material = _referenceGeometry->get(latticeR);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(&material, 1, _referenceGeometry->getLoadBalancer().rank(latticeR[0]));
#endif
  return material;
}

} // namespace opti

} // namespace olb

#endif
