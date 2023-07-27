/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;

using T = float;
using DESCRIPTOR = D3Q27<descriptors::FORCE, FreeSurface::MASS, FreeSurface::EPSILON, FreeSurface::CELL_TYPE, FreeSurface::CELL_FLAGS, FreeSurface::TEMP_MASS_EXCHANGE, FreeSurface::PREVIOUS_VELOCITY>;

/*
 * Helper since there are a lot of values to set and giving a reference to this object is easier than
 * defining the function calls accordingly
 */
struct FreeSurfaceAppHelper {
  std::array<T,3> area{{0.03, 0.03, 0.025}};
  std::array<T,3> gravity_force = {{0.,0., -9.81}};

  T char_phys_length = 0.03;
  T char_phys_vel = 0.1;
  bool has_surface_tension = true;
  T surface_tension_coefficient = 0.0661;

  std::array<T,3> initial_falling_speed{{0.,0., -6.03}};
};

template <typename T, typename DESCRIPTOR>
class FreeSurfaceFallingDrop3D final : public AnalyticalF3D<T,T> {
private:
  T lattice_size;
  // This value doesn't depend on the dimension of the descriptor. It's always 3
  std::array<T, 3> cell_values;
public:
  FreeSurfaceFallingDrop3D(T lattice_size_, const std::array<T,3>& cell_vals):AnalyticalF3D<T,T>{1}, lattice_size{lattice_size_}, cell_values{cell_vals}{}

  bool operator()(T output[], const T x[]) override {
    output[0] = cell_values[0];
    T radius = 0.00155;

    if(x[2] <= radius){
      output[0] = cell_values[2];
    }else if(x[2] <= radius + lattice_size * 1.5){
      output[0] = cell_values[1];
    }

    std::array<T,DESCRIPTOR::d> point = {0.015, 0.015, 2 * radius + lattice_size * 4};
    std::array<T,DESCRIPTOR::d> diff = {x[0] - point[0], x[1] - point[1], x[2]-point[2]};

    if((diff[0]*diff[0] + diff[1] * diff[1] + diff[2]*diff[2]) <= radius*radius){
      output[0] = cell_values[2];
    }else{
      for(int i = -1; i <= 1; ++i){
        for(int j = -1; j <= 1; ++j){
          for(int k = -1; k <= 1; ++k){
            std::array<T,DESCRIPTOR::d> shifted_diff = {diff[0]+i*lattice_size*T{1.1}, diff[1]+j*lattice_size*T{1.1}, diff[2] + k * lattice_size*T{1.1}};
            if((shifted_diff[0]*shifted_diff[0] + shifted_diff[1] * shifted_diff[1]+shifted_diff[2]*shifted_diff[2]) <= radius*radius){
              output[0] = cell_values[1];
              return true;
            }
          }
        }
      }
    }

    return true;
  }
};

template <typename T, typename DESCRIPTOR>
class FallingDropVel3D final : public AnalyticalF<3,T,T> {
private:
  T lattice_size;
  std::array<T,DESCRIPTOR::d> lattice_speed;
public:
  FallingDropVel3D(T lattice_size_, const std::array<T,DESCRIPTOR::d>& lattice_speed_):AnalyticalF<3,T,T>{3}, lattice_size{lattice_size_}, lattice_speed{lattice_speed_}{}

  bool operator()(T output[], const T x[]) override {


    T radius = 0.00155;
    std::array<T,DESCRIPTOR::d> point = {0.015, 0.015, 2 * radius + lattice_size * 2};
    std::array<T,DESCRIPTOR::d> diff = {x[0] - point[0], x[1] - point[1], x[2]-point[2]};

    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
    for(int i = -1; i <= 1; ++i){
      for(int j = -1; j <= 1; ++j){
        for(int k = -1; k <= 1; ++k){
          std::array<T,DESCRIPTOR::d> shifted_diff = {diff[0]+i*lattice_size*T{1.1}, diff[1]+j*lattice_size*T{1.1},diff[2]+k*lattice_size*T{1.1}};
          if((shifted_diff[0]*shifted_diff[0] + shifted_diff[1] * shifted_diff[1]+shifted_diff[2]*shifted_diff[2]) <= radius*radius){
            output[0] = lattice_speed[0];
            output[1] = lattice_speed[1];
            output[2] = lattice_speed[2];
            return true;
          }
        }
      }
    }

    return true;
  }
};

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1,1} );

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareFallingDrop(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry, T lattice_size, const FreeSurfaceAppHelper& helper)
{
  AnalyticalConst3D<T,T> zero( 0. );
  AnalyticalConst3D<T,T> one( 1. );
  AnalyticalConst3D<T,T> two( 2. );
  AnalyticalConst3D<T,T> four( 4. );
  FreeSurfaceFallingDrop3D<T,DESCRIPTOR> cells_analytical{ lattice_size, {0., 1., 2.}};
  FreeSurfaceFallingDrop3D<T,DESCRIPTOR> mass_analytical{ lattice_size, {0., 0.5, 1.}};

  AnalyticalConst3D<T,T> force_zero{0., 0., 0.};

  for (int i: {0,1,2}) {
    sLattice.defineField<FreeSurface::MASS>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_FLAGS>(superGeometry, i, zero);
    sLattice.defineField<descriptors::FORCE>(superGeometry, i, force_zero);
    sLattice.defineField<FreeSurface::PREVIOUS_VELOCITY>(superGeometry, i, force_zero);
  }

  sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, 1, cells_analytical);
  sLattice.defineField<FreeSurface::MASS>(superGeometry, 1, mass_analytical);
  sLattice.defineField<FreeSurface::EPSILON>(superGeometry, 1, mass_analytical);

  for (int i: {0,2}) {
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, one);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, four);
  }

  T force_factor = 1./ converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst3D<T,T> force_a{helper.gravity_force[0] * force_factor, helper.gravity_force[1] * force_factor, helper.gravity_force[2] * force_factor};
  sLattice.defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({1}), force_a);

}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry, T lattice_size, const FreeSurfaceAppHelper& helper) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>>( superGeometry, 1);
  // Material=2 -->no-slip boundary
  sLattice.defineDynamics<BounceBack<T,DESCRIPTOR>>( superGeometry, 2);
  //setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLattice.setParameter<collision::LES::Smagorinsky>(T(0.2));

  prepareFallingDrop(converter, sLattice, superGeometry, lattice_size, helper);
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                      SuperGeometry<T,3>& sGeometry,
                      T lattice_length,
                      UnitConverter<T,DESCRIPTOR> const& converter,
                      const FreeSurfaceAppHelper& helper){
  OstreamManager clout( std::cout,"setInitialValues" );

  std::array<T,DESCRIPTOR::d> lattice_speed;
  for(size_t i = 0; i < DESCRIPTOR::d; ++i){
    lattice_speed[i] = helper.initial_falling_speed[i] * converter.getPhysDeltaT() / converter.getPhysDeltaX();
  }

  FallingDropVel3D<T,DESCRIPTOR> u{lattice_length, lattice_speed};
  AnalyticalConst3D<T,T> one(1.);

  sLattice.defineRhoU( sGeometry.getMaterialIndicator({0,1,2}), one, u );
  for (int i: {0,1,2}) {
    sLattice.iniEquilibrium( sGeometry, i, one, u );
  }

  // Set up free surface communicator stages
  FreeSurface::initialize(sLattice);
  // Make the lattice ready for simulation
  sLattice.initialize();
}

namespace {
FreeSurfaceAppHelper free_surface_config;

class FreeSurfaceConfig {
public:
  T viscosity = 1e-4;
  T density = 1e3;
  T physTime = 0.01;
  T latticeRelaxationTime = .516;
  int N = 100;

  // Anti jitter value
  T transitionThreshold = 1e-3;
  // When to remove lonely cells
  T lonelyThreshold = 1.0;
};

}

void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "fallingDrop3d" );
  const int vtmIter  = converter.getLatticeTime( FreeSurfaceConfig{}.physTime /  50. );
  const int statIter = converter.getLatticeTime( FreeSurfaceConfig{}.physTime / 100. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    SuperLatticeExternalScalarField3D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon( sLattice );
    SuperLatticeExternalScalarField3D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> cells( sLattice );
    SuperLatticeExternalScalarField3D<T, DESCRIPTOR, FreeSurface::MASS> mass( sLattice );
    epsilon.getName() = "epsilon";
    cells.getName() = "cell_type";
    mass.getName() = "mass";
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.addFunctor( epsilon );
    vtmWriter.addFunctor( cells );
    vtmWriter.addFunctor( mass );

    vtmWriter.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }
}



int main(int argc, char **argv)
{
  olbInit(&argc, &argv, false, false);

  FreeSurfaceConfig c;
  OstreamManager clerr( std::cerr, "main" );
  OstreamManager clout( std::cout, "main" );

  singleton::directories().setOutputDir("./tmp/");

  FreeSurfaceAppHelper& helper = free_surface_config;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {c.N},     // resolution: number of voxels per charPhysL
    (T)   c.latticeRelaxationTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   helper.char_phys_length,     // charPhysLength: reference length of simulation geometry
    (T)   helper.char_phys_vel,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   c.viscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   c.density     // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("free surface");

  T lattice_size = helper.char_phys_length / c.N;

  T force_conversion_factor = 1./converter.getConversionFactorForce()*converter.getConversionFactorMass();

  // Convert kg / s^2
  // Basically it is multiplied with s^2 / kg = s^2 * m^3 / (kg * m^2 * m) = 1. / (velocity_factor^2 * density * length_factor)
  T surface_tension_coefficient_factor = std::pow(converter.getConversionFactorTime(),2)/ (c.density * std::pow(converter.getConversionFactorLength(),3));

  clout<<"Surface: "<<surface_tension_coefficient_factor * helper.surface_tension_coefficient<<std::endl;
  clout<<"Lattice Size: "<<converter.getPhysDeltaX()<<std::endl;

  // === 2nd Step: Prepare Geometry ===
  Vector<T,3> extend( helper.area[0], helper.area[1], helper.area[2] );
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), noOfCuboids );

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  clout<<"Overlap: "<<sLattice.getOverlap()<<std::endl;

  prepareLattice( converter, sLattice, superGeometry, lattice_size, helper);

  FreeSurface3DSetup<T,DESCRIPTOR> free_surface_setup{sLattice};

  free_surface_setup.addPostProcessor();

  // Set variables from freeSurfaceHelpers.h
  sLattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(true);
  sLattice.setParameter<FreeSurface::TRANSITION>(c.transitionThreshold);
  sLattice.setParameter<FreeSurface::LONELY_THRESHOLD>(c.lonelyThreshold);
  sLattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(helper.has_surface_tension);
  sLattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(surface_tension_coefficient_factor * helper.surface_tension_coefficient);
  sLattice.setParameter<FreeSurface::FORCE_CONVERSION_FACTOR>(force_conversion_factor);
  sLattice.setParameter<FreeSurface::LATTICE_SIZE>(converter.getPhysDeltaX());

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( c.physTime ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  setInitialValues(sLattice, superGeometry, lattice_size, converter, helper);

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( c.physTime ); ++iT ) {
    getResults( sLattice, converter, iT, superGeometry, timer );
    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();

  return 0;
}
