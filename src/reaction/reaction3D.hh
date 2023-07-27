/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 the OpenLB project
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
 * Groups all the generic 3D template files in the explicitFiniteDifference directory.
 */

#include "reactingSpecies3D.hh"
#include "rate.hh"
#include "reactionPostProcessor3D.hh"
#include "explicitFiniteDifference/explicitFiniteDifference3D.hh"

#include "eul2Lagr/eul2LagrDensity3D.hh"
#include "eul2Lagr/eul2LagrOperation3D.hh"
#include "eul2Lagr/eul2LagrPostProcessor3D.hh"
