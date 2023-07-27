/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_STAGES_H
#define CORE_STAGES_H

namespace olb {

namespace stage {

/// Communication prior to collision
struct PreCollide      { };
/// Communication after collision
struct PostCollide     { };
/// Communication after propagation
struct PostStream      { };
/// Communication after applying the post processors
struct PostPostProcess { };
/// Communication prior to coupling
struct PreCoupling     { };
/// Communication after coupling
struct PostCoupling    { };
/// On-demand communication at SuperLattice::communicate
struct Full            { };

/// Coupling post processors
struct Coupling { };

/// Collision stage
struct Collide { };

}

}

#endif
