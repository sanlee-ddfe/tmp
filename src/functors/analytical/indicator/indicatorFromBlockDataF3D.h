/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef INDICATOR_FROM_BLOCKDATA_F_3D_H
#define INDICATOR_FROM_BLOCKDATA_F_3D_H

namespace olb {

template <typename S>
class IndicatorBlockData3D : public IndicatorF3D<S> {
private:
  BlockData<3,S,S> _blockData;
  S _deltaR;
  bool _invert;
public:
  /// constructor
  IndicatorBlockData3D(BlockData<3,S,S> blockData, Vector<S,3> extend, Vector<S,3> origin, S deltaR, bool invert=false);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]) override;
};


}

#endif

