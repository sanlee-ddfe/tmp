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
 * This file contains functions used for the particle-wall and particle-particle communication.
*/

#ifndef PARTICLE_CONTACT_COMMUNICATION_FUNCTIONS_H
#define PARTICLE_CONTACT_COMMUNICATION_FUNCTIONS_H

namespace olb {

namespace particles {

namespace contact {

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void communicateContacts(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "contactCommunication");

  if (singleton::mpi().getSize() > 1) {
    int sumParticleContact = 0;
    int sumWallContact     = 0;
    int countParticleContactsPerRank[singleton::mpi().getSize()];
    int countWallContactsPerRank[singleton::mpi().getSize()];
    int localCountParticleContacts = 0;
    int localCountWallContacts     = 0;

    //clout << "communicate number of contacts ..." << std::endl;
    if (singleton::mpi().isMainProcessor()) {
      countParticleContactsPerRank[singleton::mpi().bossId()] =
          contactContainer.particleContacts.size();
      sumParticleContact +=
          countParticleContactsPerRank[singleton::mpi().bossId()];
      countWallContactsPerRank[singleton::mpi().bossId()] =
          contactContainer.wallContacts.size();
      sumWallContact += countWallContactsPerRank[singleton::mpi().bossId()];
      for (int rank = 0; rank < singleton::mpi().getSize(); ++rank) {
        if (rank != singleton::mpi().bossId()) {
          singleton::mpi().receive(&countParticleContactsPerRank[rank], 1,
                                   rank);
          sumParticleContact += countParticleContactsPerRank[rank];
          singleton::mpi().receive(&countWallContactsPerRank[rank], 1, rank);
          sumWallContact += countWallContactsPerRank[rank];
        }
      }
    }
    else {
      localCountParticleContacts =
          contactContainer.particleContacts
              .size(); // size of std::vector (called contacts in container)
      singleton::mpi().send(&localCountParticleContacts, 1,
                            singleton::mpi().bossId()); // 0 is rank of master
      localCountWallContacts = contactContainer.wallContacts.size();
      singleton::mpi().send(&localCountWallContacts, 1,
                            singleton::mpi().bossId());
    }
    //clout << "communicate number of contacts ... OK" << std::endl;

    //clout << "broadcasting number of contacts ..." << std::endl;
    singleton::mpi().bCast(&sumParticleContact, 1, singleton::mpi().bossId());
    singleton::mpi().bCast(&sumWallContact, 1, singleton::mpi().bossId());
    //std::cout << sumParticleContact << std::endl;
    //singleton::mpi().barrier();
    //clout << "broadcasting number of contacts ... OK" << std::endl;

    if (sumParticleContact > 0 || sumWallContact > 0) {
      ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> container(
          sumParticleContact, sumWallContact);
      int particleContactObjSize = sizeof(PARTICLECONTACTTYPE);
      int wallContactObjSize     = sizeof(WALLCONTACTTYPE);

      //clout << "Bundling contact objects ..." << std::endl;
      if (singleton::mpi().isMainProcessor()) {
        PARTICLECONTACTTYPE* currParticle = container.particleContacts.data();
        WALLCONTACTTYPE*     currWall     = container.wallContacts.data();
        for (int rank = 0; rank < singleton::mpi().getSize(); ++rank) {
          std::uint8_t* bytesParticle = (std::uint8_t*)currParticle;
          std::uint8_t* bytesWall     = (std::uint8_t*)currWall;
          if (rank != singleton::mpi().bossId()) {
            singleton::mpi().receive(bytesParticle,
                                     countParticleContactsPerRank[rank] *
                                         particleContactObjSize,
                                     rank);
            singleton::mpi().receive(
                bytesWall, countWallContactsPerRank[rank] * wallContactObjSize,
                rank);
          }
          else {
            std::memcpy(
                bytesParticle,
                (std::uint8_t*)contactContainer.particleContacts.data(),
                countParticleContactsPerRank[singleton::mpi().bossId()] *
                    particleContactObjSize);
            std::memcpy(bytesWall,
                        (std::uint8_t*)contactContainer.wallContacts.data(),
                        countWallContactsPerRank[singleton::mpi().bossId()] *
                            wallContactObjSize);
          }
          currParticle += countParticleContactsPerRank[rank];
          currWall += countWallContactsPerRank[rank];
        }
      }
      else {
        PARTICLECONTACTTYPE* currParticle =
            contactContainer.particleContacts.data();
        WALLCONTACTTYPE* currWall      = contactContainer.wallContacts.data();
        std::uint8_t*    bytesParticle = (std::uint8_t*)currParticle;
        std::uint8_t*    bytesWall     = (std::uint8_t*)currWall;
        singleton::mpi().send(
            bytesParticle, localCountParticleContacts * particleContactObjSize,
            singleton::mpi().bossId());
        singleton::mpi().send(bytesWall,
                              localCountWallContacts * wallContactObjSize,
                              singleton::mpi().bossId());
      }
      //clout << "Bundling contact objects ... OK" << std::endl;

      if (singleton::mpi().isMainProcessor()) {
        container.combineContacts();
      }

      //clout << "Broadcasting contact objects ..." << std::endl;
      std::uint8_t* containerParticlesBytes =
          (std::uint8_t*)container.particleContacts.data();
      singleton::mpi().bCast(containerParticlesBytes,
                             sumParticleContact * particleContactObjSize,
                             singleton::mpi().bossId());
      std::uint8_t* containerWallBytes =
          (std::uint8_t*)container.wallContacts.data();
      singleton::mpi().bCast(containerWallBytes,
                             sumWallContact * wallContactObjSize,
                             singleton::mpi().bossId());
      //singleton::mpi().barrier();
      //clout << "Broadcasting contact objects ... OK" << std::endl;
      /*
            clout << "#######" << std::endl;
            clout << "id: " << container.wallContacts[3].id << std::endl;
            clout << "#######" << std::endl;
      */
      //std::cout << "size before: " << container.particleContacts.size() << std::endl;
      container.cleanContacts();
      //std::cout << "size after: " << container.particleContacts.size() << std::endl;
      contactContainer = container;
    }
  }
#endif
}

} //namespace contact

} //namespace particles

} //namespace olb

#endif
