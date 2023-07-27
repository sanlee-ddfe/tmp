/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Thomas Henn, Mathias J. Krause
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

#ifndef MATERIALSTLBOUNDARY3D_H
#define MATERIALSTLBOUNDARY3D_H

#include <set>
#include <vector>

#include "geometry/superGeometry.h"
#include "particles/subgrid3DLegacyFramework/particleSystem3D.h"
#include "boundary3D.h"

namespace olb {



template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;


/*
 * Particle boundary based on the fluids material number using signed distance function.
 * If a particle moves in the vincinity of a lattice node with
 * specified material number it is set inactive and its velocity is set to 0.
 **/

template<typename T, template<typename U> class PARTICLETYPE>
class MaterialSTLBoundary3D : public Boundary3D<T, PARTICLETYPE> {
public:
  /// Constructor
  MaterialSTLBoundary3D(SuperGeometry<T,3>& sg);
  /// Copy constructor
  MaterialSTLBoundary3D(MaterialSTLBoundary3D<T, PARTICLETYPE>& f);
  /// Constructor with set of material numbers
  /*MaterialBoundary3D(SuperGeometry<T,3>& sg,
                     std::set<int> material);*/
  MaterialSTLBoundary3D(SuperGeometry<T,3>& sg,
                     std::set<int> material, STLreader<T> &stlReader);
/*MaterialBoundary3D(SuperGeometry<T,3>& sg,
                     std::set<int> material);*/
  ~MaterialSTLBoundary3D() override
  {
  }
  /// Add a single material number
  void addMaterial(int mat)
  {
    _materials.insert(mat);
  }
  /// Add several material numbers
  void addMaterial(std::vector<int> mats)
  {
    for (unsigned i=0; i< mats.size(); ++i) {
      _materials.insert(mats[i]);
    }
  }
  /// Apply the boundary condition
  void applyBoundary(
    typename std::deque<PARTICLETYPE<T> >::iterator& p,
    ParticleSystem3D<T, PARTICLETYPE>& psSys) override;

private:
  SuperGeometry<T,3>& _sg;
  std::set<int> _materials;
  std::set<int>::iterator _matIter;
  STLreader<T> &_stlReader;
  //  int x0, y0, z0;
};

}

#endif /* MATERIALSTLBOUNDARY3D_H */
