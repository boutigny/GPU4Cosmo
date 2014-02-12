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

# CUDA code generation flags
# 
GENCODE_SM10  := -gencode arch=compute_10,code=sm_10
GENCODE_SM20  := -gencode arch=compute_20,code=sm_20
GENCODE_SM30  := -gencode arch=compute_30,code=sm_30 \
                 -gencode arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_FLAGS := $(GENCODE_SM10) $(GENCODE_SM20) $(GENCODE_SM30)

# Source files
#
.SUFFIXES : .o .cu .cpp

CUDA_SOURCES := angCorr.cu
CUDA_OBJECTS  = $(CUDA_SOURCES:.cu=.o)

CXX_SOURCES  := cppRef.cpp 
CXX_OBJECTS   = $(CXX_SOURCES:.cpp=.o)

# Targets
#
all: ../bin/angCorr ../bin/cppRef

../bin/angCorr: $(CUDA_OBJECTS) 
	@echo -e "\nLinking" $@ 
	$(NVCC) $(ALL_LD_FLAGS) $(CUDA_OBJECTS) -o $@

../bin/cppRef: $(CXX_OBJECTS)
	@echo -e "\nLinking" $@
	$(CXX) $(ROOT_LD_FLAGS) $(CXX_OBJECTS) -o $@

.cu.o:
	@echo -e "\nCompiling" $<
	$(NVCC) $(ALL_CC_FLAGS) $(GENCODE_FLAGS) -o $@ -c $<

.cpp.o:
	@echo -e "\nCompiling" $<
	$(CXX) $(CXX_FLAGS) $(ROOT_CXX_FLAGS) -o $@ -c $<

clean:
	-rm -f angCorr cppRef *.o *~

buildall: clean all
