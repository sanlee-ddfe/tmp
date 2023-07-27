/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2021 Mathias Krause
 *                2022      Adrian Kummerlaender
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

#ifndef BLOCK_LATTICE_HH
#define BLOCK_LATTICE_HH

#include "blockLattice.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
BlockLattice<T,DESCRIPTOR>::BlockLattice(Vector<int,DESCRIPTOR::d> size, int padding, Platform platform)
  : BlockStructure<DESCRIPTOR>(size, padding),
    _platform(platform),
    _statisticsEnabled{true},
    _statistics{nullptr}
{ }

template<typename T, typename DESCRIPTOR>
BlockLattice<T,DESCRIPTOR>::~BlockLattice()
{
  if (_statistics) {
    delete _statistics;
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::initialize()
{
  addPostProcessor(typeid(stage::PostStream), meta::id<StatisticsPostProcessor>());

  // Todo: Move to StatisticsPostProcessor::setup
  _statistics = new LatticeStatistics<T>;
  _statistics->initialize();

  setProcessingContext(ProcessingContext::Simulation);
  postProcess();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineRho(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        rho(&rhoTmp,physR);
        get(loc).defineRho(rhoTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineU(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp,physR);
        get(loc).defineU(uTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineRhoU(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
  AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        rho(&rhoTmp,physR);
        u(uTmp,physR);
        get(loc).defineRhoU(rhoTmp,uTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& popF)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T pop[DESCRIPTOR::q];
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        popF(pop,physR);
        get(loc).definePopulations(pop);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::definePopulations(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                   BlockF<T,DESCRIPTOR::d>& popF)
{
  if (!indicator.isEmpty()) {
    T pop[DESCRIPTOR::q];
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        popF(pop, loc);
        get(loc).definePopulations(pop);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  if (!indicator.isEmpty()) {
    // Don't use FieldD here as long as AnalyticalF is fixed to T
    std::array<T,DESCRIPTOR::template size<FIELD>()> fieldTmp;
    T physR[DESCRIPTOR::d] = { };
    this->forSpatialLocationsParallel([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        field(fieldTmp.data(), physR);
        get(loc).template setField<FIELD>(fieldTmp.data());
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                             BlockF<T,DESCRIPTOR::d>& field)
{
  if (!indicator.isEmpty()) {
    FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        field(fieldTmp.data(), loc.data());
        get(loc).template setField<FIELD>(fieldTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicatorF,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  BlockIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicator(indicatorF, blockGeometry);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        get(loc).iniEquilibrium(rhoTmp, uTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                BlockF<T,DESCRIPTOR::d>& u)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, loc.data());
        rho(&rhoTmp, physR);
        get(loc).iniEquilibrium(rhoTmp, uTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniRegularized(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    T piTmp[util::TensorVal<DESCRIPTOR>::n] = { };
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        pi(piTmp, physR);
        get(loc).iniRegularized(rhoTmp, uTmp, piTmp);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineDynamics(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, Dynamics<T,DESCRIPTOR>* dynamics)
{
  if (!indicator.isEmpty()) {
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        defineDynamics(loc, dynamics);
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void BlockLattice<T,DESCRIPTOR>::defineDynamics(BlockIndicatorF<T,DESCRIPTOR::d>& indicator)
{
  if (!indicator.isEmpty()) {
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        setDynamics(this->getCellId(loc),
                    DynamicsPromise(meta::id<DYNAMICS>{}));
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics)
{
  this->forCellIndices([&](CellID iCell) {
    setDynamics(iCell, dynamics);
  });
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::collide()
{
  if (_customCollisionO) {
    _customCollisionO->operator()(*this);
  } else {
    if constexpr (isPlatformCPU(PLATFORM)) {
      _dynamicsMap.collide(CollisionDispatchStrategy::Dominant);
    } else {
      _dynamicsMap.collide(CollisionDispatchStrategy::Individual);
    }
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::stream()
{
  if constexpr (PLATFORM == Platform::GPU_CUDA) {
    DESCRIPTOR::template filter<descriptors::is_propagatable_field>
              ::for_each([&](auto field) {
      auto& population = getField(field.get());
      for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        population[iPop].rotate(this->getNeighborDistance(descriptors::c<DESCRIPTOR>(iPop)));
      }
      // Required for resolving populations in (gather,scatter)_any_fields kernel
      getDataRegistry().refreshDeviceFieldArray(population);
    });
  } else {
    DESCRIPTOR::template filter<descriptors::is_propagatable_field>
              ::for_each([&](auto field) {
      auto& population = getField(field.get());
      for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        population[iPop].rotate(this->getNeighborDistance(descriptors::c<DESCRIPTOR>(iPop)));
      }
    });
  }
}

/// Operator for striping off density offset
/**
 * Used by BlockLattice::stripeOffDensityOffset
 **/
struct StripeOffDensityOffsetO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct OFFSET : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<OFFSET>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto offset = parameters.template get<OFFSET>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] -= descriptors::t<V,DESCRIPTOR>(iPop) * offset;
    }
  }
};

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  if (!hasPostProcessor(typeid(StripeOffDensityOffsetO),
                        meta::id<StripeOffDensityOffsetO>{})) {
    addPostProcessor(typeid(StripeOffDensityOffsetO),
                     meta::id<StripeOffDensityOffsetO>{});
  }
  // Update offset and execute
  getData<OperatorParameters<StripeOffDensityOffsetO>>()
    .template set<StripeOffDensityOffsetO::OFFSET>(offset);
  postProcess(typeid(StripeOffDensityOffsetO));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::addLatticeCoupling(LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners)
{
  addPostProcessor<stage::Coupling>(lcGen.generate(partners));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::addLatticeCoupling(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                    LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners)
{
  if (!indicator.isEmpty()) {
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (this->getNeighborhoodRadius(loc) >= 1) {
        if (indicator(loc)) {
          std::unique_ptr<LatticeCouplingGenerator<T,DESCRIPTOR>> extractedLcGen{ lcGen.clone() };
          if (extractedLcGen->extract(0, 0)) {
            extractedLcGen->shift(loc);
            addLatticeCoupling(*extractedLcGen, partners);
          }
        }
      }
    });
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::executeCoupling()
{
  postProcess<stage::Coupling>();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLattice<T,DESCRIPTOR>::getStatistics()
{
  return *_statistics;
}

template<typename T, typename DESCRIPTOR>
const LatticeStatistics<T>& BlockLattice<T,DESCRIPTOR>::getStatistics() const
{
  return *_statistics;
}


template<typename T, typename DESCRIPTOR, Platform PLATFORM>
ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::ConcreteBlockLattice(Vector<int,DESCRIPTOR::d> size, int padding)
  : BlockLattice<T,DESCRIPTOR>(size, padding, PLATFORM),
    _data(),
    _descriptorFields(),
    _dynamicsMap(*this)
{
  DESCRIPTOR::fields_t::template for_each([&](auto id) {
    using field = typename decltype(id)::type;
    using field_type = Array<field>;
    auto& fieldArray = _data.template allocate<field_type>(this->getNcells());
    _descriptorFields.template set<field>(&fieldArray);
    _communicatables[typeid(field)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename field::template column_type<T>,PLATFORM>::type,
                 DESCRIPTOR::template size<field>()>
    >(fieldArray));
    _data.template setSerialization<field_type>(true);
  });

  Dynamics<T,DESCRIPTOR>* noDynamics = _dynamicsMap.get(DynamicsPromise(meta::id<NoDynamics<T,DESCRIPTOR>>{}));
  for (CellID iCell=0; iCell < this->getNcells(); ++iCell) {
    _dynamicsMap.set(iCell, noDynamics);
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
bool ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::hasData() const
{
  return _data.template provides<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
const auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getData() const
{
  OLB_ASSERT(_data.template provides<FIELD_TYPE>(),
             "FIELD_TYPE must be allocated to be accessed");
  return _data.template get<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getData()
{
  if (_data.template provides<FIELD_TYPE>()) {
    return _data.template get<FIELD_TYPE>();
  } else {
    // TODO: Implement more generic approach to constructing arbitrary data from specific args
    auto& data = _data.template allocate<FIELD_TYPE>(this->getNcells());
    // Manage serializables and communicatables for array fields
    using concrete_data_t = typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>;
    if constexpr (std::is_base_of_v<ColumnVectorBase,concrete_data_t>) {
      using field_t = typename concrete_data_t::field_t;
      if constexpr (field_t::isSerializable()) {
        _data.template setSerialization<FIELD_TYPE>(true);
      }
      _communicatables[typeid(field_t)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
        ColumnVector<typename ImplementationOf<typename field_t::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<field_t>()>
      >(data));
    }
    return data;
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getField(FIELD)
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
const auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getField(FIELD) const
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<const FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::addPostProcessor(
  std::type_index stage, PostProcessor<T,DESCRIPTOR>* postProcessor)
{
  auto [postProcessorsOfPriority, _] = _postProcessors[stage].try_emplace(postProcessor->getPriority(), this);
  std::get<1>(*postProcessorsOfPriority).addLegacy(postProcessor);
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::addPostProcessor(
  std::type_index stage, const PostProcessorGenerator<T,DESCRIPTOR>& ppGen)
{
  addPostProcessor(stage, ppGen.generate());
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::addPostProcessor(
  std::type_index stage, BlockIndicatorF<T,DESCRIPTOR::d>& indicator, const PostProcessorGenerator<T,DESCRIPTOR>& ppGen)
{
  if (!indicator.isEmpty()) {
    this->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        std::unique_ptr<PostProcessorGenerator<T,DESCRIPTOR>> extractedPpGen{ ppGen.clone() };
        if (extractedPpGen->extract(0, 0)) {
          extractedPpGen->shift(loc);
          addPostProcessor(stage, *extractedPpGen);
        }
      }
    });
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::postProcess(std::type_index stage)
{
  for (auto& [_, postProcessorsOfPriority] : _postProcessors[stage]) {
    postProcessorsOfPriority.apply();
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  return _data.getNblock();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  return _data.getSerializableSize();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data, loadingMode);

  return dataPtr;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::postLoad()
{
  auto& population = getField<descriptors::POPULATION>();
  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    population[iPop].postLoad();
  }
}


}

#endif
