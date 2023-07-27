/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
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
 * Groups all the include files for 2D dynamics in
 * the dataStructure directory.
 */

#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "loadBalancer.h"
#include "blockLoadBalancer.h"
#include "heuristicLoadBalancer.h"
#include "randomLoadBalancer.h"
#include "mpiManager.h"
#include "mpiManagerAD.hh"  // includes aDiff, but “S” is not defined -> it is working if you use the right order in the main cf. apps/mathias/bifurcation-pi
#include "ompManager.h"
#include "superStructure.h"
#include "blockCommunicator.h"
#include "superCommunicator.h"
#include "blockCommunicationNeighborhood.h"
#include "superCommunicationTagCoordinator.hh"
#include "mpiGroup.h"

#endif
