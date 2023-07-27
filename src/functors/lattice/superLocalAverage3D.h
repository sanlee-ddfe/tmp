/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_LOCAL_AVERAGE_3D_H
#define SUPER_LOCAL_AVERAGE_3D_H

#include "superBaseF3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "utilities/functorPtr.h"

namespace olb {


/// Averages given functor inside the local sphere
template <typename T, typename W = T>
class SuperLocalAverage3D final : public SuperF3D<T,W> {
private:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
  const T _radius;

public:
  /// Primary constructor
  /**
   * \param f          Functor to be locally averaged
   * \param indicatorF Indicator describing relevant cells
   * \param radius     Radius of the locality sphere
   **/
  SuperLocalAverage3D(FunctorPtr<SuperF3D<T>>&&          f,
                      FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                      T radius);

  /**
   * Returns average of functor \e _f evaluated on all cells both inside a sphere
   * of \e _radius around \e input and indicated by \e _indicatorF.
   **/
  bool operator() (W output[], const int input[]) override;
};


}

#endif
