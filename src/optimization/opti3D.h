/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Mathias J. Krause, Benjamin Förster
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
 * Groups all the 3D include files for the directory opti.
 */
#include "utilities/aDiff.h"
#include "adjointLbSolver.h"
#include "controlledFunctions3D.h"
#include "dualDynamics.h"
#include "dualFunctors3D.h"
#include "dualMrtDynamics.h"
#include "optiCase.h"
#include "optiCaseAD.h"
#include "optiCaseDual.h"
#include "optimizer.h"
#include "optimizerBarzilaiBorwein.h"
#include "optimizerLineSearch.h"
#include "optimizerLBFGS.h"
#include "optimizerSteepestDecent.h"
#include "optiSolverParameters.h"
#include "dualLatticeDescriptors.h"
#include "sequentialInterpolationF3D.h"
