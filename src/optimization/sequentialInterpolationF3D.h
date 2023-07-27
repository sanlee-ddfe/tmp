/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster, Adrian Kummerlaender
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

/* \file
 * In order to get the unparallelized functors from controlledFunctions3D.h
 * running, the old interpolation functors have been copied from 8915a03ff.
 * TODO: parallelize the functors in opti code s.t. this file becomes
 * unnecassary.
 */

#ifndef SEQUENTIAL_INTERPOLATION_F_3D_H
#define SEQUENTIAL_INTERPOLATION_F_3D_H

#include "functors/analytical/analyticalF.h"
#include "functors/lattice/blockBaseF3D.h"
#include "functors/lattice/superBaseF3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/blockGeometry.h"
#include "geometry/superGeometry.h"

namespace olb {


/// Converts block functors to analytical functors (special)
template <typename T, typename W = T>
class SpecialSequentialAnalyticalFfromBlockF3D final : public AnalyticalF3D<T,W> {
protected:
  BlockF3D<W>& _f;
  Cuboid3D<T>& _cuboid;
  Vector<T,3> _delta;
  T _scale;
public:
  SpecialSequentialAnalyticalFfromBlockF3D(BlockF3D<W>& f, Cuboid3D<T>& cuboid, Vector<T,3> delta, T scale = 1.);
  bool operator() (W output[], const T physC[]) override;
};

/// Converts block functors to analytical functors
template <typename T, typename W = T>
class SequentialAnalyticalFfromBlockF3D final : public AnalyticalF3D<T,W> {
protected:
  BlockF3D<W>& _f;
  Cuboid3D<T>& _cuboid;
  const int    _overlap;
public:
  SequentialAnalyticalFfromBlockF3D(BlockF3D<W>& f, Cuboid3D<T>& cuboid, const int overlap);
  bool operator() (W output[], const T physC[]) override;
};

/// Converts super functors to analytical functors
template <typename T, typename W = T>
class SequentialAnalyticalFfromSuperF3D final : public AnalyticalF3D<T,W> {
protected:
  const bool _communicateToAll;
  const bool _communicateOverlap;

  SuperF3D<T,W>&       _f;
  CuboidGeometry3D<T>& _cuboidGeometry;
  int                  _overlap;

  std::vector<std::unique_ptr<SequentialAnalyticalFfromBlockF3D<T,W>>> _blockF;
public:
  SequentialAnalyticalFfromSuperF3D(SuperF3D<T,W>& f,
                          bool communicateToAll=false,
                          int overlap=-1,
                          bool communicateOverlap=true);
  bool operator() (W output[], const T physC[]) override;

  /// \return Size of _blockF vector
  int getBlockFSize() const;
  /// \return _blockF[iCloc]
  SequentialAnalyticalFfromBlockF3D<T,W>& getBlockF(int iCloc);
};


} // end namespace olb

#endif
