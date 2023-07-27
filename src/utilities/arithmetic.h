/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef UTILITY_ARITHMETIC_H
#define UTILITY_ARITHMETIC_H

#include "utilities/omath.h"
#include <functional>

namespace olb {

namespace util {

/// Wrapper of function object std::minus with special handling for bool
/**
 * \tparam T Domain of the substraction operation, _without_ operation is
 *           performed for boolean inputs.
 *
 * Note that specialization is not required for boolean union (add) and
 * intersection (multiply) functors as their behavior is implicitly defined
 * by the corresponding arithmetic operation on the integer representation.
 **/
template <typename T>
struct minus {
  /// symbol character for functor naming
  static const char symbol = '-';

  constexpr T operator() (const T& lhs, const T& rhs) const
  {
    return std::minus<T>()(lhs, rhs);
  }
};

/// Operator specialization for boolean _without_ operation
/**
 * i.e. implements what one reasonably expects to happen when substracting
 *      indicator functors.
 **/
template <>
constexpr bool minus<bool>::operator() (const bool& lhs, const bool& rhs) const
{
  return lhs && !rhs;
}

/// Wrapper of function object std::plus
template <typename T>
struct plus : public std::plus<T> {
  /// symbol character for functor naming
  static const char symbol = '+';
};

/// Wrapper of function object std::multiplies
template <typename T>
struct multiplies : public std::multiplies<T> {
  /// symbol character for functor naming
  static const char symbol = '*';
};

/// Wrapper of function object std::divides
template <typename T>
struct divides : public std::divides<T> {
  /// symbol character for functor naming
  static const char symbol = '/';
};

/// Power function object
template <typename T>
struct power {
  /// symbol character for functor naming
  static const char symbol = '^';

  constexpr T operator() (const T& base, const T& exponent) const
  {
    return util::pow(base, exponent);
  }
};

/// Wrapper function object for util::max
// imitates stl definition of std::plus
template<class T = void> struct maxOp {
  /// symbol character for functor naming
  static const char symbol = 'M';

  constexpr T operator()(const T& x, const T& y) const {
    return util::max(x, y);
  }
};

template<> struct maxOp<void> {
  /// symbol character for functor naming
  static const char symbol = 'M';

  template<class T, class U> constexpr auto operator()(T&& t, U&& u) const
    -> decltype(util::max(std::forward<T>(t), std::forward<U>(u))) {
      return util::max(std::forward<T>(t), std::forward<U>(u));
    }

  using is_transparent = void;
};

/// Wrapper function object for util::min
// imitates stl definition of std::plus
template<class T = void> struct minOp {
  /// symbol character for functor naming
  static const char symbol = 'm';

  constexpr T operator()(const T& x, const T& y) const {
    return util::min(x, y);
  }
};

template<> struct minOp<void> {
  /// symbol character for functor naming
  static const char symbol = 'm';

  template<class T, class U> constexpr auto operator()(T&& t, U&& u) const
    -> decltype(util::min(std::forward<T>(t), std::forward<U>(u))) {
      return util::min(std::forward<T>(t), std::forward<U>(u));
    }

  using is_transparent = void;
};

} // end namespace util

} // end namespace olb

#endif
