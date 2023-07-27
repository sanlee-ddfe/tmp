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

#ifndef PARTICLE_DYNAMICS_BASE_H
#define PARTICLE_DYNAMICS_BASE_H

namespace olb {

namespace particles {

namespace dynamics {


/// Basic particle dynamics
template<typename T, typename PARTICLETYPE>
struct ParticleDynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~ParticleDynamics() { }
  /// Implementation of the processing step
  virtual void process(Particle<T,PARTICLETYPE>& particle, T timeStepSize) =0;
  /// read and write access to name
  std::string& getName();
  /// read only access to name
  std::string const& getName() const;
private:
  std::string _name;
};

/// No particle dynamics equivalent to no lattice dynamics
template<typename T, typename PARTICLETYPE>
class NoParticleDynamics : public ParticleDynamics<T,PARTICLETYPE> {
public:
  NoParticleDynamics( T rhoDummy );
  /// Processing step
  void process(Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Standard dynamics for particles
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::active_particles>
class VerletParticleDynamics : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamics( );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Verlet particle dynamics only considering translation (ignoring rotation)
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsTranslationOnly : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamicsTranslationOnly( );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Verlet particle dynamics only considering rotation (ignoring translation)
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsRotationOnly : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamicsRotationOnly( );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Standard dynamics with wall reflection
template<typename T, typename PARTICLETYPE, bool useCubicBounds=false,
         typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsVelocityWallReflection
  : public VerletParticleDynamics<T,PARTICLETYPE,PCONDITION> {
public:
  /// Constructor
  VerletParticleDynamicsVelocityWallReflection( SolidBoundary<T,PARTICLETYPE::d>& solidBoundary );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  SolidBoundary<T,PARTICLETYPE::d>& _solidBoundary;
};

/// Standard dynamics with wall capture
template<typename T, typename PARTICLETYPE, bool useCubicBounds=false,
         typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsWallCapture
  : public VerletParticleDynamics<T,PARTICLETYPE,PCONDITION> {
public:
  /// Constructor
  VerletParticleDynamicsWallCapture( SolidBoundary<T,PARTICLETYPE::d>& solidBoundary );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  SolidBoundary<T,PARTICLETYPE::d>& _solidBoundary;
};

/// Standard dynamics with material capture
/// - implies necessity of SuperGeometry, so no DEM only possible
template<typename T, typename PARTICLETYPE,
         typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsMaterialCapture
  : public VerletParticleDynamics<T,PARTICLETYPE,PCONDITION> {
public:
  /// Constructor
  VerletParticleDynamicsMaterialCapture(
    SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& _materialIndicator;
};

/// Standard dynamics with wall capture and material number checks
/// - implies necessity of SuperGeometry, so no DEM only possible
template<typename T, typename PARTICLETYPE,
         typename PCONDITION=conditions::active_particles>
class VerletParticleDynamicsMaterialAwareWallCapture
  : public VerletParticleDynamics<T,PARTICLETYPE,PCONDITION> {
public:
  /// Constructor
  VerletParticleDynamicsMaterialAwareWallCapture( SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
    SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  SolidBoundary<T,PARTICLETYPE::d>& _solidBoundary;
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& _materialIndicator;
};


/// Verlet dynamics for particles aware of their DYNAMIC_STATE
/// - For now, manly used for wall-flow application
///   - source: https://gitlab.com/openlb/olb/-/blob/private/nicolas/main4/apps/nicolas/A_WallFlow/wallFlow/README.md
///   - partially assuming axis aligned surfaces and provision of main flow direction
template<typename T, typename PARTICLETYPE>
class ParticleDetachmentDynamics
  : public VerletParticleDynamics<T,PARTICLETYPE,conditions::active_particles> {
public:
  /// Constructor
  ParticleDetachmentDynamics( SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
                              Vector<T,PARTICLETYPE::d>& mainFlowDirection,
                              T tiltThreshold = 0.3*M_PI );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  SolidBoundary<T,PARTICLETYPE::d>& _solidBoundary;
  Vector<T,PARTICLETYPE::d> _mainFlowDirection;
  T _tiltThreshold; //Tilt angle initiating particle release
};




//TODO: REWORK FOLLOWING DYNAMICS ACCORDING TO OBOVE ONES


//TODO: remove for release
/// Velocity verlet particle dynamics with limitation of position and velocity by checking domain bounds
/// in cartesion direcion and simple adhesive force threshold allowing particles only to move
/// when both a normal and tangential force threshold have been surpassed
template<typename T, typename PARTICLETYPE>
class VerletParticleDynamicsCubicBoundsAdhesion : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamicsCubicBoundsAdhesion( PhysR<T,PARTICLETYPE::d>& domainMin,
      PhysR<T,PARTICLETYPE::d>& domainMax );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  PhysR<T,PARTICLETYPE::d> _domainMin;
  PhysR<T,PARTICLETYPE::d> _domainMax;
};


//TODO: remove for release
/// Velocity verlet particle dynamics with deposition modelling by checking domain bounds
/// in cartesion direcion
template<typename T, typename PARTICLETYPE, typename DEPOSITION_MODEL>
class VerletParticleDynamicsCubicBoundsDeposition : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamicsCubicBoundsDeposition( PhysR<T,PARTICLETYPE::d>& domainMin,
    PhysR<T,PARTICLETYPE::d>& domainMax, DEPOSITION_MODEL& depositionModel );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  PhysR<T,PARTICLETYPE::d> _domainMin;
  PhysR<T,PARTICLETYPE::d> _domainMax;
  DEPOSITION_MODEL& _depositionModel;
};

} //namespace dynamics

} //namespace particles

} //namespace olb

#endif
