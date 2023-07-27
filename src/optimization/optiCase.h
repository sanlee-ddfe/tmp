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


#ifndef OPTI_CASE_H
#define OPTI_CASE_H

#include <functional>

#include "io/xmlReader.h"

namespace olb {

namespace opti {


/// Abstract base class for optimization tasks
// provides function evaluation and gradient computation
template <typename S, typename C>
class OptiCase{
public:
  virtual S evaluateObjective(const C& control, unsigned optiStep=0) = 0;
  virtual void computeDerivatives (
    const C& control, C& derivatives, unsigned optiStep=0) = 0;
};


/// Gradient is just passed as a function (and not computed by an own routine)
template <typename S, typename C>
class OptiCaseAnalytical : public OptiCase<S,C> {

protected:
  std::function<S (const C&)>         _function;
  std::function<void (const C&, C&)>  _derivative;

public:
  OptiCaseAnalytical(std::function<S (const C&)> function,
    std::function<void (const C&, C&)> derivative)
   : _function(function), _derivative(derivative)
  { }

  S evaluateObjective(
    const C& control, unsigned optiStep=0) override {
    return _function(control);
  }

  void computeDerivatives(
    const C& control, C& derivatives, unsigned optiStep=0) override {
    _derivative(control, derivatives);
  }
};

template <typename S, typename C>
OptiCaseAnalytical(std::function<S (const C&)>,
                   std::function<void (const C&, C&)>)
  -> OptiCaseAnalytical<S,C>;


// Abstract base class for gradient computation with difference quotients
template <typename S, typename C>
class OptiCaseDQ : public OptiCase<S,C> {

protected:
  std::function<S (const C&)>          _function;
  S                                    _stepWidth {1.e-8};

  bool                                 _objectiveComputed {false};
  S                                    _objective;

public:
  explicit OptiCaseDQ(std::function<S (const C&)> function) :
    _function(function)
  { }

  OptiCaseDQ(std::function<S (const C&)> function, S stepWidth) :
    _function(function), _stepWidth(stepWidth)
  { }

  S evaluateObjective(
    const C& control, unsigned optiStep=0) override {
    _objective = _function(control);
    _objectiveComputed = true;
    return _objective;
  }
};


/// Gradient computation with forward difference quotients
template <typename S, typename C>
class OptiCaseFDQ : public OptiCaseDQ<S,C> {
private:
  OstreamManager                                  clout {std::cout, "OptiCaseFDQ"};

public:
  explicit OptiCaseFDQ(std::function<S (const C&)> function) :
    OptiCaseDQ<S,C>(function)
  { }

  OptiCaseFDQ(std::function<S (const C&)> function, S stepWidth) :
    OptiCaseDQ<S,C>(function, stepWidth)
  { }

  void computeDerivatives(const C& control,
    C& derivatives, unsigned optiStep=0) override
  {
    assert((control.size() == derivatives.size()));

    if (!(this->_objectiveComputed)) {
      this->evaluateObjective(control);
    }
    const S objective(this->_objective);

    for (std::size_t it = 0; it < control.size(); ++it)
    {
      C shiftedControl(control);
      shiftedControl[it] += this->_stepWidth;
      S shiftedObjective = this->evaluateObjective(shiftedControl);
      derivatives[it] = (shiftedObjective - objective) / this->_stepWidth;
    }
    this->_objectiveComputed = false;
  }
};

template <typename S, typename C>
OptiCaseFDQ(std::function<S (const C&)>)
  -> OptiCaseFDQ<S,C>;

template <typename S, typename C>
OptiCaseFDQ(std::function<S (const C&)>, S)
  -> OptiCaseFDQ<S,C>;



/// Gradient computation with central difference quotients
template <typename S, typename C>
class OptiCaseCDQ : public OptiCaseDQ<S,C> {
private:
  OstreamManager                                  clout {std::cout, "OptiCaseCDQ"};

public:
  explicit OptiCaseCDQ(std::function<S (const C&)> function) :
    OptiCaseDQ<S,C>(function)
  { }

  OptiCaseCDQ(std::function<S (const C&)> function, S stepWidth) :
    OptiCaseDQ<S,C>(function, stepWidth)
  { }

  void computeDerivatives(const C& control,
    C& derivatives, unsigned optiStep=0) override
  {
    assert((control.size() == derivatives.size()));

    for (std::size_t it = 0; it < control.size(); ++it)
    {
      C shiftedControl(control);
      shiftedControl[it] += this->_stepWidth;
      S shiftedObjective_plus = this->evaluateObjective(shiftedControl);

      shiftedControl[it] = control[it] - this->_stepWidth;
      S shiftedObjective_minus = this->evaluateObjective(shiftedControl);

      derivatives[it] = 0.5 * (shiftedObjective_plus - shiftedObjective_minus) / this->_stepWidth;
    }
  }
};

template <typename S, typename C>
OptiCaseCDQ(std::function<S (const C&)>)
  -> OptiCaseCDQ<S,C>;

template <typename S, typename C>
OptiCaseCDQ(std::function<S (const C&)>, S)
  -> OptiCaseCDQ<S,C>;


} // namespace opti

} // namespace olb

#endif
