
NVCC=nvcc
CXX=g++

ROOTLIB=/usr/lib64/root/
ROOTINC=/usr/include/root/
HEALPIXLIB=/usr/lib64/
HEALPIXINC=/usr/include/healpix/

LIBDIR=-L${ROOTLIB} -L${HEALPIXLIB} 
LIBS=-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -m64 -lhealpix_cxx -lhealpix_cxxsupport -lhealpix_fft -lcfitsio 
INCLUDES=-I${ROOTINC} -I${HEALPIXINC} 

ROOTFLAGS=$(shell root-config --libs --cflags)

all: bin/angCorr bin/angCorr.nogpu

clean:
	rm -f bin/angCorr bin/angCorr.nogpu

bin/angCorr: 2ptCorr/angCorr.cu
	$(NVCC) 2ptCorr/angCorr.cu $(LIBDIR) $(LIBS) $(INCLUDES) -o bin/angCorr

bin/angCorr.nogpu: 2ptCorr/cppRef.cpp
	$(CXX) 2ptCorr/cppRef.cpp $(ROOTFLAGS) $(LIBDIR) $(LIBS) $(INCLUDES) -o bin/angCorr.nogpu
