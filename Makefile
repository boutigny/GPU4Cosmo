
NVCC=nvcc

ROOTLIB=/usr/lib64/root/
ROOTINC=/usr/include/root/
HEALPIXLIB=/usr/lib64/
HEALPIXINC=/usr/include/healpix/

LIBDIR=-L${ROOTLIB} -L${HEALPIXLIB} 
LIBS=-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -m64 -lhealpix_cxx -lhealpix_cxxsupport -lhealpix_fft -lcfitsio 
INCLUDES=-I${ROOTINC} -I${HEALPIXINC} 

all: bin/angCorr

clean:
	rm -f bin/angCorr

bin/angCorr: 2ptCorr/angCorr.cu
	$(NVCC) 2ptCorr/angCorr.cu $(LIBDIR) $(LIBS) $(INCLUDES) -o bin/angCorr