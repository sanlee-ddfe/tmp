# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Tested using CUDA 11.4 and OpenMPI 4.1 (CUDA aware)
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Adjust CUDA_ARCH to match your specifc GPU
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun -np 2 ./cavity3d` (All processes share default GPU, not optimal)
#
# Usage on a multi GPU system: (recommended when using MPI, use non-MPI version on single GPU systems)
#  - Run `mpirun -np 4 bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./cavity3d'
#    (for a 4 GPU system, further process mapping advisable, consult cluster documentation)
#
# (CUDA_)CXXFLAGS and (CUDA_)LDFLAGS may need to be adjusted depending on the specific system
# environment.
#
# One advantage of mixed compilation is that `make no-cuda-recompile` is available,
# significantly speeding up the compilation if the set of used dynamics / post
# processors has not been expanded (i.e. if no new GPU kernels need costly
# instantiation).

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++17

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

# Compiler to use for CUDA-enabled files
CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++17
# Adjust to enable resolution of libcuda, libcudart, libcudadevrt
CUDA_LDFLAGS    := -L/run/opengl-driver/lib
# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
