/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012, 2021 Mathias J. Krause, Julius Jeßberger
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
* An OptiCase using Algorithmic Differentiation (AD)
*/


#ifndef OPTI_CASE_AD_H
#define OPTI_CASE_AD_H


#include "utilities/aDiff.h"
#include "optiCase.h"
#include "utilities/vectorHelpers.h"


namespace olb {

namespace opti {

template<typename S, std::size_t n,
  template<typename> typename C = util::StdVector>
class OptiCaseAD : public OptiCase<S,C<S>> {

protected:
  using T = util::ADf<S,n>;
  std::function<S (const C<S>&)>      _function;
  std::function<T (const C<T>&)>      _adFunction;

public:
  explicit OptiCaseAD() = default;

  OptiCaseAD(
    std::function<S (const C<S>&)> function,
    std::function<T (const C<T>&)> adFunction) :
    _function(function), _adFunction(adFunction)
  { }

  S evaluateObjective(
    const C<S>& control, unsigned optiStep=0) override
  {
    return _function(control);
  }

  void computeDerivatives(const C<S>& control,
    C<S>& derivatives, unsigned optiStep=0) override
  {
    assert((control.size() == derivatives.size()));

    auto adControl = util::iniAD<n,S,C>(control);

    T adResult = _adFunction(adControl);

    for(std::size_t it = 0; it < control.size(); ++it){
      derivatives[it] = adResult.d(it);
    }
  }
};


/** Interface for OptiCaseAD that performs Lattice-Boltzmann-Solver
 * construction itself (from xml file)
 */
template<typename S, std::size_t n, template<typename> typename SOLVER,
  template<typename> typename C = util::StdVector>
class OptiCaseAdForSolver : public OptiCaseAD<S,n,C> {

protected:
  using typename OptiCaseAD<S,n,C>::T;
  std::shared_ptr<SOLVER<S>>    _solver;
  std::shared_ptr<SOLVER<T>>    _adSolver;

public:
  OptiCaseAdForSolver(XMLreader const& xml) :
    _solver(createLbSolver<SOLVER<S>> (xml)),
    _adSolver(createLbSolver<SOLVER<T>> (xml))
  {
    this->_function = getCallable<S>(_solver);
    this->_adFunction = getCallable<T>(_adSolver);
  }
};

}

}

#endif
