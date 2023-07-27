/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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

/** \file A fluid flow simulation test case - execution script
 * Cf. testFlowOpti3d.h for implementation details
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "testFlowOpti3d.h"

using namespace olb;
using namespace olb::parameters;

using S = FLOATING_POINT_TYPE;
constexpr unsigned numberOfDerivatives (3);
using T = util::ADf<S,numberOfDerivatives>;

XMLreader config("parameter.xml");


/** Compute sensitivities with forward difference quotients
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void sensitivitiesFDQ() {
  OstreamManager clout(std::cout, "sensitivitiesFDQ");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  testFlow->parameters(Output()).verbose = false;
  OptiCaseFDQ optiCase(getCallable<S>(testFlow));

  // evaluate function
  std::vector<S> control (numberOfDerivatives, 2.5);
  clout << "objective = " << optiCase.evaluateObjective(control) << std::endl;

  // compute derivatives
  std::vector<S> derivatives (numberOfDerivatives, 12.);
  optiCase.computeDerivatives(control, derivatives);
  util::print(derivatives, "derivatives", clout, ' ');
}

/** Compute sensitivities with automatic differentiation
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void sensitivitiesAD() {
  // preparation
  OstreamManager clout(std::cout, "sensitivitiesAD");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  testFlow->parameters(Output()).verbose = false;
  testFlowAD->parameters(Output()).verbose = false;
  testFlowAD->parameters(VisualizationVTK()).output = false;

  OptiCaseAD<S, numberOfDerivatives> optiCase(
    getCallable<S>(testFlow),
    getCallable<T>(testFlowAD));

  // evaluate function
  std::vector<S> control (numberOfDerivatives, 2.5);
  util::print(control, "control", clout, ' ');
  clout << "objective = " << optiCase.evaluateObjective(control) << std::endl;

  // compute (forward AD) derivatives
  std::vector<S> derivatives (numberOfDerivatives, 12.);
  optiCase.computeDerivatives(control, derivatives);
  util::print(derivatives, "derivatives", clout, ' ');
}

/** Variant for evaluation + forward AD sensitivity computation with arbitrary
 * number of control variables: "control" scales the force cuboid-wise
 * (independent of the cuboid decomposition which is used for parallelization).
 */
void sensitivitiesAD_variant(unsigned numberOfControls)
{
  // preparation
  OstreamManager clout(std::cout, "sensitivitiesAD");

  std::vector<S> control (numberOfControls, 2.5);
  util::print(control, "control", clout, ' ');

  // evaluate function
  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  auto function = getCallable<S> (testFlow);
  clout << "objective = " << function(control) << std::endl;

  // compute (forward AD) derivatives
  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  auto derivative = getCallable<T> (testFlowAD);

  auto controlAD = util::copyAs<T,S,util::StdVector>(control);
  util::iniDiagonal(&controlAD[0], numberOfDerivatives);
  auto resultAD = derivative(controlAD);
  util::print(resultAD.d(), "derivatives", clout, ' ');
}

/** Optimization with automatic differentiation
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void optiAD() {
  OstreamManager clout(std::cout, "optiAD");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);

  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  testFlow->parameters(Output()).verbose = false;
  testFlowAD->parameters(Output()).verbose = false;
  testFlowAD->parameters(VisualizationVTK()).output = false;

  bool discreteObjective = testFlow->parameters(Opti()).optiReferenceMode;
  auto referenceSolver = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  std::vector<S> referenceControl (numberOfDerivatives, 1.0);
  if (discreteObjective) {
    // compute numerical solution for velocity/ dissipation for
    // correct control variables
    referenceSolver->parameters(Opti()).computeObjective = false;
    referenceSolver->parameters(Output()).verbose = false;
    referenceSolver->parameters(VisualizationVTK()).output = false;

    referenceSolver->parameters(Opti()).applyControl(referenceControl);
    referenceSolver->solve();

    testFlow->parameters(Opti()).referenceSolution
     = referenceSolver->parameters(Results()).solution;

    testFlowAD->parameters(Opti()).referenceSolution
     = referenceSolver->parameters(Results()).solution;
  }

  OptiCaseAD<S, numberOfDerivatives> optiCase(
    getCallable<S>(testFlow),
    getCallable<T>(testFlowAD));

  auto optimizer = createOptimizerLBFGS<S>(config, numberOfDerivatives);

  std::vector<S> startValues (numberOfDerivatives, 0.0);
  optimizer->setControl(startValues);
  optimizer->setReferenceControl(referenceControl);

  // Execute optimization
  optimizer->optimize(optiCase);

  clout << "objective = " << optimizer->getObjective() << std::endl;
  util::print(optimizer->getControl(),    "control   ", clout, ' ');
  util::print(optimizer->getDerivative(), "derivative", clout, ' ');
}

/** Optimization with adjoint LBM
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void optiAdjoint() {
  OstreamManager clout(std::cout, "optiAdjoint");

  OptiCaseDual<S,TestFlowSolverOptiAdjoint> optiCase (config);
  auto optimizer = createOptimizerLBFGS<S>(config, optiCase._dimCtrl);

  S startValue;
  config.readOrWarn<S>("Optimization", "StartValue", "", startValue);
  optimizer->setStartValue(optiCase.getInitialControl(startValue));

  optimizer->setReferenceControl(optiCase.getReferenceControl());

  optimizer->optimize(optiCase);
}


int main(int argc, char **argv)
{
  olbInit(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  clout << "\nCompute sensitivities with forward difference quotients" << std::endl;
  sensitivitiesFDQ();

  clout << "\nCompute sensitivities with AD" << std::endl;
  sensitivitiesAD();

  clout << "\nCompute sensitivities with AD for 5 blocks" << std::endl;
  sensitivitiesAD_variant(5);

  clout << "\nOptimize with AD" << std::endl;
  optiAD();

  clout << "\nOptimize with adjoint lbm" << std::endl;
  optiAdjoint();
}
