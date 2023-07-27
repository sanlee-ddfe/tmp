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

#ifndef CPU_SISD_MASK_H_
#define CPU_SISD_MASK_H_

#include "column.h"

namespace olb {

template <typename T>
class ConcreteBlockMask<T,Platform::CPU_SISD> final : public Serializable {
private:
  cpu::sisd::Column<bool> _mask;
  std::size_t _weight;

public:
  ConcreteBlockMask(std::size_t size):
    _mask(size),
    _weight(0)
  { }

  bool operator[](std::size_t i) const {
    return _mask[i];
  }

  void set(std::size_t i, bool active) {
           if ( _mask[i] && !active) {
      _weight -= 1;
    } else if (!_mask[i] &&  active) {
      _weight += 1;
    }
    _mask[i] = active;
  }

  std::size_t weight() const {
    return _weight;
  }

  void setProcessingContext(ProcessingContext) { }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T>
std::size_t ConcreteBlockMask<T,Platform::CPU_SISD>::getNblock() const
{
  return 1 + _mask.getNblock();
}

template<typename T>
std::size_t ConcreteBlockMask<T,Platform::CPU_SISD>::getSerializableSize() const
{
  return sizeof(_weight) + _mask.getSerializableSize();
}

template<typename T>
bool* ConcreteBlockMask<T,Platform::CPU_SISD>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _weight);
  registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _mask, loadingMode);

  return dataPtr;
}


}

#endif
