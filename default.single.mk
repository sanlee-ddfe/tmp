include $(OLB_ROOT)/rules.mk

all: dependencies core $(EXAMPLE)

INCLUDE_DIRS := $(OLB_ROOT)/src

ifeq ($(USE_EMBEDDED_DEPENDENCIES), ON)
INCLUDE_DIRS += \
	$(OLB_ROOT)/external/zlib \
	$(OLB_ROOT)/external/tinyxml

LDFLAGS := -L$(OLB_ROOT)/external/lib $(LDFLAGS)

dependencies:
	$(MAKE) -C $(OLB_ROOT)/external
else
.PHONY: dependencies
endif

LDFLAGS += -L$(OLB_ROOT)/build/lib

core:
	$(MAKE) -C $(OLB_ROOT) core

CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
DEP_FILES := $(CPP_FILES:.cpp=.d)

INCLUDE_FLAGS := -I$(subst $(EMPTY) $(EMPTY), -I,$(INCLUDE_DIRS))

%.d: $(CPP_FILES)
	@$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

$(EXAMPLE): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $@ -lolbcore $(LDFLAGS)

.PHONY: onlysample
onlysample: $(EXAMPLE)

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONE: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(EXAMPLE)

-include $(DEP_FILES)
