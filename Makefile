###########################################################################
## OpenLB Makefile

DEFAULT_TARGETS :=
CLEAN_TARGETS :=

# Call DEFAULT_TARGETS for plain `make`
all: default

###########################################################################
## Default build configuration

# Include config.mk environment (optional)
-include config.mk

include rules.mk

###########################################################################
## Embedded dependencies (optional)

ifeq ($(USE_EMBEDDED_DEPENDENCIES), ON)
dependencies:
	$(MAKE) CXX=$(CXX) CC=$(CC) -C external

clean-dependencies:
	$(MAKE) CXX=$(CXX) CC=$(CC) -C external clean

CLEAN_TARGETS += clean-dependencies
else
.PHONY: dependencies
endif

DEFAULT_TARGETS += dependencies

###########################################################################
## Core library

CORE_CPP_FILES := \
	src/communication/mpiManager.cpp \
	src/communication/ompManager.cpp \
	src/core/olbInit.cpp \
	src/io/ostreamManager.cpp
CORE_OBJ_FILES := $(CORE_CPP_FILES:.cpp=.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fPIC -Isrc/ -c $< -o $@

build/lib:
	mkdir -p build/lib

build/lib/libolbcore.a: build/lib $(CORE_OBJ_FILES)
	ar rc $@ $(CORE_OBJ_FILES)

core: build/lib/libolbcore.a

DEFAULT_TARGETS += core

clean-core:
	rm -f $(CORE_OBJ_FILES)
	rm -f build/lib/libolbcore.a

CLEAN_TARGETS += clean-core

###########################################################################
## Examples

EXAMPLES := $(dir $(shell find examples -name 'Makefile'))

CLEAN_SAMPLES_TARGETS :=

$(EXAMPLES): dependencies core
	$(MAKE) -C $@ onlysample

.PHONY: $(EXAMPLES)

samples: $(EXAMPLES)

define clean_directory
clean-$(1):
	+$(MAKE) -C $(1) clean

CLEAN_SAMPLES_TARGETS += clean-$(1)
endef

$(foreach sample,$(EXAMPLES), \
  $(eval $(call clean_directory,$(sample))))

clean-samples: $(CLEAN_SAMPLES_TARGETS)

###########################################################################
## Code generation (CSE)

CSE_GENERATORS := $(shell find src/ -type f -name '*.cse.h.template')
CSE_GENERATEES := $(patsubst %.cse.h.template, %.cse.h, $(CSE_GENERATORS))

%.cse.h: %.cse.h.template
	$(eval $@_GUARD := $(shell grep -Po "#ifndef\ \K[A-Z_]*_CSE_H" $<))
	python codegen/cse.py $< $($@_GUARD) $@

cse: $(CSE_GENERATEES)

###########################################################################
## Doxygen documentation

doxygen:
	doxygen doc/DoxygenConfig

###########################################################################
## Cleaning

clean: $(CLEAN_TARGETS)

clean-all: clean-samples clean

###########################################################################
## Default targets (called for plain `make`)

default: $(DEFAULT_TARGETS)
