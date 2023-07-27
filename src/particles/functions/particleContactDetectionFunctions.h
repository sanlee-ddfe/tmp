/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt, Mathias J. Krause
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

/*
 * This file contains functions used for the particle-wall and particle-particle contact detection.
*/

#ifndef PARTICLE_CONTACT_DETECTION_FUNCTIONS_H
#define PARTICLE_CONTACT_DETECTION_FUNCTIONS_H

namespace olb {

namespace particles {

namespace contact {

template <typename CONTACTTYPE>
auto getContactIterator(std::vector<CONTACTTYPE>&                     contacts,
                        const std::function<bool(const CONTACTTYPE&)> condition)
{
  auto contactIt = std::find_if(contacts.begin(), contacts.end(),
                                [&condition](const auto& contact) -> bool {
                                  return condition(contact);
                                });
  return contactIt;
}

std::array<std::size_t, 2>
sortParticleIDs(const std::array<std::size_t, 2>& ids)
{
  return std::array<std::size_t, 2>(
      {util::min(ids[0], ids[1]), util::max(ids[0], ids[1])});
}

template <typename T, typename PARTICLETYPE, typename F>
void unifyPositions(Particle<T, PARTICLETYPE>&       particle1,
                    Particle<T, PARTICLETYPE>&       particle2,
                    const PhysR<T, PARTICLETYPE::d>& cellMin,
                    const PhysR<T, PARTICLETYPE::d>& cellMax,
                    F                                getSetupPeriodicity,
                    std::array<PhysR<T, PARTICLETYPE::d>, 2>& positions,
                    T                                         deltaX)
{
  using namespace descriptors;

  const PhysR<T, PARTICLETYPE::d> pos1 =
      particles::access::getPosition(particle1);
  const PhysR<T, PARTICLETYPE::d> pos2 =
      particles::access::getPosition(particle2);
  positions[0] = pos1;
  positions[1] = pos2;
}

template <typename T, typename PARTICLETYPE, typename F>
PhysR<T, PARTICLETYPE::d>
unifyPosition(Particle<T, PARTICLETYPE>&       particle,
              const PhysR<T, PARTICLETYPE::d>& cellMin,
              const PhysR<T, PARTICLETYPE::d>& cellMax, F getSetupPeriodicity,
              T deltaX)
{
  using namespace descriptors;
  const PhysR<T, PARTICLETYPE::d> pos =
      particles::access::getPosition(particle);
  PhysR<T, PARTICLETYPE::d> unifiedPosition(pos);

  return unifiedPosition;
}

template <typename T, typename PARTICLETYPE, typename F>
PhysR<T, PARTICLETYPE::d>
evalContactPosition(Particle<T, PARTICLETYPE>&       particle,
                    const PhysR<T, PARTICLETYPE::d>& particlePos,
                    const PhysR<T, PARTICLETYPE::d>& contactPos,
                    const PhysR<T, PARTICLETYPE::d>& cellMin,
                    const PhysR<T, PARTICLETYPE::d>& cellMax,
                    F getSetupPeriodicity, T deltaX)
{
  return contactPos;
}

template <typename PARTICLECONTACTTYPE, bool IS_INPUT_SORTED = false>
bool particleContactConsistsOfIDs(PARTICLECONTACTTYPE& particleContact,
                                    const std::array<size_t, 2>& ids)
{
  if constexpr (!IS_INPUT_SORTED) {
    return (particleContact.getIDs()[0] == ids[0] &&
            particleContact.getIDs()[1] == ids[1]) ||
           (particleContact.getIDs()[0] == ids[1] &&
            particleContact.getIDs()[1] == ids[0]);
  }
  else {
    return particleContact.getIDs()[0] == ids[0] &&
           particleContact.getIDs()[1] == ids[1];
  }
  __builtin_unreachable();
}

template <typename T, typename PARTICLETYPE, bool CONVEX, typename F>
void updateContact(ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d,
                                                             CONVEX>& contact,
                   Particle<T, PARTICLETYPE>&                         particle1,
                   Particle<T, PARTICLETYPE>&                         particle2,
                   const PhysR<T, PARTICLETYPE::d>& contactPos,
                   const PhysR<T, PARTICLETYPE::d>& cellMin,
                   const PhysR<T, PARTICLETYPE::d>& cellMax,
                   F getSetupPeriodicity, T deltaX)
{
  if (!contact.isParticlePositionUpdated()) {
    std::array<PhysR<T, PARTICLETYPE::d>, 2> positions;
    unifyPositions(particle1, particle2, cellMin, cellMax, getSetupPeriodicity,
                   positions, deltaX);
    contact.setParticlePositions(positions);
  }
  // It doesn't matter which particle is used, since the contact must be in both
  contact.updateMinMax(evalContactPosition(
      particle1, contact.getParticlePosition(contact.getIDs()[0]),
      contactPos, cellMin, cellMax, getSetupPeriodicity, deltaX));
}

template <typename T, typename PARTICLETYPE, bool CONVEX, typename F>
void updateContact(WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d,
                                                         CONVEX>& contact,
                   Particle<T, PARTICLETYPE>&                     particle,
                   const PhysR<T, PARTICLETYPE::d>&               contactPos,
                   const PhysR<T, PARTICLETYPE::d>&               cellMin,
                   const PhysR<T, PARTICLETYPE::d>&               cellMax,
                   F getSetupPeriodicity, T deltaX)
{
  if (!contact.isParticlePositionUpdated()) {
    contact.setParticlePosition(
        unifyPosition(particle, cellMin, cellMax, getSetupPeriodicity, deltaX));
  }
  contact.updateMinMax(
      evalContactPosition(particle, contact.getParticlePosition(), contactPos,
                          cellMin, cellMax, getSetupPeriodicity, deltaX));
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, typename F>
void updateContacts(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    //const std::array<size_t, 2>& ids,
    std::array<std::size_t, 2>&& ids, const PhysR<T, PARTICLETYPE::d>& pos,
    Particle<T, PARTICLETYPE>& particle1, Particle<T, PARTICLETYPE>& particle2,
    const PhysR<T, PARTICLETYPE::d>& cellMin,
    const PhysR<T, PARTICLETYPE::d>& cellMax, F getSetupPeriodicity, T deltaX)
{
  // find the contact by checking particle ids
  const std::function<bool(const PARTICLECONTACTTYPE&)> condition =
      [&ids](const PARTICLECONTACTTYPE& contact) {
        // we assume that ids aren't sorted
        return particleContactConsistsOfIDs(contact, ids);
      };
  auto contactIt =
      getContactIterator(contactContainer.particleContacts, condition);

  if (contactIt != std::end(contactContainer.particleContacts)) {
    if (contactIt->isNew()) {
      if (particleContactConsistsOfIDs<PARTICLECONTACTTYPE, true>(
              *contactIt, ids)) {
        updateContact(*contactIt, particle1, particle2, pos, cellMin, cellMax,
                      getSetupPeriodicity, deltaX);
      }
      else {
        updateContact(*contactIt, particle2, particle1, pos, cellMin, cellMax,
                      getSetupPeriodicity, deltaX);
      }
    }
  }
  else {
    contactContainer.particleContacts.push_back(PARTICLECONTACTTYPE(ids));
    if (particleContactConsistsOfIDs<PARTICLECONTACTTYPE, true>(
            *(contactContainer.particleContacts.end() - 1), ids)) {
      updateContact(*(contactContainer.particleContacts.end() - 1), particle1,
                    particle2, pos, cellMin, cellMax, getSetupPeriodicity,
                    deltaX);
    }
    else {
      updateContact(*(contactContainer.particleContacts.end() - 1), particle2,
                    particle1, pos, cellMin, cellMax, getSetupPeriodicity,
                    deltaX);
    }
  }
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, typename F>
void updateContacts(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    size_t particleID, unsigned wallID, const PhysR<T, PARTICLETYPE::d>& pos,
    Particle<T, PARTICLETYPE>&       particle,
    const PhysR<T, PARTICLETYPE::d>& cellMin,
    const PhysR<T, PARTICLETYPE::d>& cellMax, F getSetupPeriodicity, T deltaX)
{
  // find the contact by checking ids
  const std::function<bool(const WALLCONTACTTYPE&)> condition =
      [&particleID, &wallID](const WALLCONTACTTYPE& contact) {
        return particleID == contact.getParticleID() &&
               wallID == contact.getWallID();
      };
  auto contactIt =
      getContactIterator(contactContainer.wallContacts, condition);

  if (contactIt != std::end(contactContainer.wallContacts)) {
    if (contactIt->isNew()) {
      updateContact(*contactIt, particle, pos, cellMin, cellMax,
                    getSetupPeriodicity, deltaX);
    }
  }
  else {
    contactContainer.wallContacts.push_back(
        WALLCONTACTTYPE(particleID, wallID));
    updateContact(*(contactContainer.wallContacts.end() - 1), particle, pos,
                  cellMin, cellMax, getSetupPeriodicity, deltaX);
  }
}

} //namespace contact

} //namespace particles

} //namespace olb

#endif
