###########################################################################
## conditional settings

# Select std::thread foundation for use in thread pool / async job system
ifdef THREAD_POOL
ifeq ($(THREAD_POOL), PTHREAD)
	CXXFLAGS += -pthread
	LDFLAGS += -lpthread
endif
# Default to pthread
else
	CXXFLAGS += -pthread
	LDFLAGS += -lpthread
endif

ifeq ($(PARALLEL_MODE), MPI)
	CXXFLAGS += -DPARALLEL_MODE_MPI $(MPIFLAGS)
	PARALLEL_FLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS)
	LDFLAGS  := $(MPIFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
	CXXFLAGS += -DPARALLEL_MODE_OMP $(OMPFLAGS)
	PARALLEL_FLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS)
	LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
	CXXFLAGS += -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS)
	PARALLEL_FLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS)
	LDFLAGS  := $(MPIFLAGS) $(OMPFLAGS) $(LDFLAGS)
endif

LDFLAGS += $(EXTRA_LDFLAGS)

FILTERED_FEATURES := $(FEATURES)

ifneq ($(filter DISABLE_CSE,$(FILTERED_FEATURES)),)
	CXXFLAGS += -DDISABLE_CSE
	FILTERED_FEATURES := $(filter-out DISABLE_CSE,$(FILTERED_FEATURES))
endif

LDFLAGS += $(if $(filter $(FEATURES), OPENBLAS),-lopenblas)
LDFLAGS += -lz -ltinyxml

ifneq ($(filter CPU_SIMD,$(PLATFORMS)),)
	LDFLAGS += -lrt
endif

ifneq ($(filter GPU_CUDA,$(PLATFORMS)),)
## | CUDA Architecture | Version    |
## |-------------------+------------|
## | Fermi             | 20         |
## | Kepler            | 30, 35, 37 |
## | Maxwell           | 50, 52, 53 |
## | Pascal            | 60, 61, 62 |
## | Volta             | 70, 72     |
## | Turing            | 75         |
## | Ampere            | 80, 86, 87 |
	CUDA_ARCH ?= 75

	CUDA_LDFLAGS += -lcuda -lcudadevrt -lcudart

	ifndef CUDA_CXXFLAGS
		CUDA_CXXFLAGS := -O3 -std=c++17
	endif

	CUDA_CXXFLAGS += --generate-code=arch=compute_$(CUDA_ARCH),code=[compute_$(CUDA_ARCH),sm_$(CUDA_ARCH)]
	CUDA_CXXFLAGS += --extended-lambda --expt-relaxed-constexpr
	# Enable relocatable device code (only needed for mixed mode but speeds up non-mixed compilation significantly)
	CUDA_CXXFLAGS += -rdc=true
	# Hide some distracting warnings on calling device host functions from host
	CUDA_CXXFLAGS += -Xcudafe "--diag_suppress=implicit_return_from_non_void_function --display_error_number --diag_suppress=20014 --diag_suppress=20011"

	ifndef CUDA_CXX
		CXXFLAGS += --forward-unknown-to-host-compiler -x cu
		CXXFLAGS += $(CUDA_CXXFLAGS)
		LDFLAGS += $(CUDA_LDFLAGS)
	endif
endif

PLATFORM_FLAGS := $(foreach platform,$(PLATFORMS),-DPLATFORM_$(platform))
FEATURE_FLAGS := $(foreach feature,$(FILTERED_FEATURES),-DFEATURE_$(feature))

CXXFLAGS += $(PLATFORM_FLAGS)
CXXFLAGS += $(FEATURE_FLAGS)

ifdef FLOATING_POINT_TYPE
	CXXFLAGS += -DDEFAULT_FLOATING_POINT_TYPE=$(FLOATING_POINT_TYPE)
endif
