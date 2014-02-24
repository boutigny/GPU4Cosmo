# C++ compiler settings
#
CXX       ?= /usr/bin/g++
CXX_FLAGS := -Wall -Wno-unknown-pragmas

# ROOT settings
#
ROOTSYS ?= /usr/local/root/pro
ROOT_CXX_FLAGS := -I$(shell $(ROOTSYS)/bin/root-config --incdir) -m64
ROOT_LD_FLAGS  := -L$(shell $(ROOTSYS)/bin/root-config --libdir) -lCore -lCint -lRIO -lNet   \
                  -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
                  -lMathCore -lThread -lm -ldl


# HEALPIX settings
#
HEALPIX_PATH      ?= /sps/lsst/Library/Healpix_3.11
HEALPIX_CXX_FLAGS := -I$(HEALPIX_PATH)/src/cxx/basic_gcc/include
HEALPIX_LD_FLAGS  := -L$(HEALPIX_PATH)/src/cxx/basic_gcc/lib -lhealpix_cxx -lcxxsupport

# CUDA compiler settings
#
CUDA_PATH     ?= /opt/cuda-5.5
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(CXX)
NVCC_CC_FLAGS := -g -G -I$(CUDA_PATH)/include
NVCC_LD_FLAGS := -L$(CUDA_PATH)/lib -lcuda -lm -lstdc++

ALL_CC_FLAGS  := $(NVCC_CC_FLAGS)
ALL_CC_FLAGS  += $(addprefix -Xcompiler ,$(CXX_FLAGS))
ALL_CC_FLAGS  += $(addprefix -Xcompiler ,$(ROOT_CXX_FLAGS))
ALL_CC_FLAGS  += $(addprefix -Xcompiler ,$(HEALPIX_CXX_FLAGS))

ALL_LD_FLAGS  := $(addprefix -Xlinker ,$(LD_FLAGS))
ALL_LD_FLAGS  += $(addprefix -Xlinker ,$(ROOT_LD_FLAGS))
ALL_LD_FLAGS  += $(addprefix -Xlinker ,$(HEALPIX_LD_FLAGS))

# Source files
#
.SUFFIXES : .o .cu .cpp

CUDA_SOURCES := angCorr.cu cudaRef.cu
CUDA_OBJECTS  = $(CUDA_SOURCES:.cu=.o)

CXX_SOURCES  := cpuRef.cpp
CXX_OBJECTS   = $(CXX_SOURCES:.cpp=.o)

# Directories
#
BIN_DIR = ../bin
OUTPUT_DIR = ../output

# Targets
#
.PHONY: all cpuref cudaref angcorr clean buildall testcpu testcuda testall

all: cpuref cudaref angcorr

cpuref: $(BIN_DIR)/cpuref

cudaref: $(BIN_DIR)/cudaref

angcorr: $(BIN_DIR)/angcorr

$(BIN_DIR)/angcorr: angCorr.o
	@echo -e "\nLinking" $@
	$(NVCC) $(ALL_LD_FLAGS) angCorr.o -o $@

$(BIN_DIR)/cpuref: $(CXX_OBJECTS)
	@echo -e "\nLinking" $@
	$(CXX) $(ROOT_LD_FLAGS) $(CXX_OBJECTS) -o $@

$(BIN_DIR)/cudaref: cudaRef.o
	@echo -e "\nLinking" $@
	$(NVCC) $(ALL_LD_FLAGS) cudaRef.o -o $@

.cu.o:
	@echo -e "\nCompiling" $<
	$(NVCC) $(ALL_CC_FLAGS) -o $@ -c $<

.cpp.o:
	@echo -e "\nCompiling" $<
	$(CXX) $(CXX_FLAGS) $(ROOT_CXX_FLAGS) -o $@ -c $<

clean:
	-rm -f $(BIN_DIR)/angCorr $(BIN_DIR)/cpuRef $(BIN_DIR)/cudaRef *.o *~ $(OUTPUT_DIR)/*.root

buildall: clean all


testall: testcpu testcuda

testcpu:
	@echo -e "\nTesting CPU reference program" $<
	@./testCPU.sh

testcuda:
	@echo -e "\nTesting CUDA reference program" $<
	@./testCUDA.sh
