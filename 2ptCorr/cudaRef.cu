#define bcc_cxx

#include <iostream>
#include <assert.h>
#include "math.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH2F.h"

#include "bcc.h"

using namespace std;

__device__ float angDist(float ra_1, float dec_1, float ra_2, float dec_2) {
//
// Return angular distance between two points in spherical coordinate 
//
	
	float sindec1 = sinf(dec_1);
    float cosdec1 = cosf(dec_1);
    float sindec2 = sinf(dec_2);
    float cosdec2 = cosf(dec_2);
    float cosra2_ra1 = cosf(ra_2-ra_1);
    float sinra2_ra1 = sinf(ra_2-ra_1);

    float aux = (cosdec1 * sindec2) - (sindec1 * cosdec2 * cosra2_ra1);
    float num = (cosdec2 * cosdec2 * sinra2_ra1 * sinra2_ra1) + (aux * aux);
    float den = (sindec1 * sindec2) + (cosdec1 * cosdec2 * cosra2_ra1);

    return atan2f(sqrtf(num),den);
	
/*	float theta;
    theta = atan(sqrt(cos(dec_2)*cos(dec_2)*sin(ra_2-ra_1)*sin(ra_2-ra_1) + 
                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))*
                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))) /
                                (sin(dec_1)*sin(dec_2) + cos(dec_1)*cos(dec_2)*cos(ra_2-ra_1)));

    return theta; */
}

__global__ void etBim(const float *ra, const float *dec, int *histo, const int nbins, const float bin_width, const float ang_min, 
                    const float ang_max, const int index, const int n) {

    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;
    if(i < n) {
        float angle = angDist(ra[index], dec[index], ra[i+index+1], dec[i+index+1]);
        if (angle < ang_max) {
            int bin = (int)((angle-ang_min)/bin_width);
            int ibin = threadIdx.x*nbins+bin;
            histo[ibin]++;
        }
    }
}

class corelGPU {
	public:
		corelGPU(int, int, float, float);
		~corelGPU();
		void loadValues(float*, float*);
		void run();
		int* getHisto();
		size_t getBlockSize();
	private:
		int _nentries;
		int _nbins;
		float _vMin;
		float _vMax;
		static const size_t _blockSize = 1024;
		size_t _numBytes;
		size_t _numBytesHisto;
		float _binWidth;
		
		// _devX and _devY are 2 floating point arrays on device memory containing the values of the variables to be corelated
		// There is 1 array per device (2 on ccgpu)
		float* _devX[2];
		float* _devY[2];
		
		// _devHisto is an histogram containing the result of the correlation computation (distribution of the correlation values)
		// 1 array per device
		int* _devHisto[2];
};

corelGPU::corelGPU(int nentries, int nbins, float vMin, float vMax) : _nentries(nentries), _nbins(nbins),
			_vMin(vMin), _vMax(vMax) 
{
	_numBytes = nentries*sizeof(float);
	_numBytesHisto = nbins*_blockSize*sizeof(int);
	
	_devX[0] = 0;
	_devX[1] = 0;
	_devY[0] = 0;
	_devY[1] = 0;
	_devHisto[0] =0;
	_devHisto[1] = 0;
	
	_binWidth = (vMax-vMin)/(float)nbins;
    
	// allocate device memory on both GPU
    //
	for (int dev=0; dev<2; dev++) {
		cudaSetDevice(dev);
		cudaMalloc((void**)&_devX[dev], _numBytes);
		cudaMalloc((void**)&_devY[dev], _numBytes);
		cudaMalloc((void**)&_devHisto[dev], _numBytesHisto);
		cudaMemset(_devHisto[dev],0,_numBytesHisto);
	
		// check allocation
		//
		assert(_devX[dev] != 0 && _devY[dev] != 0 && _devHisto[dev] != 0);
	}
}
corelGPU::~corelGPU() {

	for (int dev=0; dev<2; dev++) {
		cudaSetDevice(dev);
		cudaFree(_devX[dev]);
		cudaFree(_devY[dev]);
		cudaFree(_devHisto[dev]);
	}
}
void corelGPU::loadValues(float* xVal, float* yVal) {

    // Copy xVal and yVal arrays to device memory
    //
	for (int dev=0; dev<2; dev++) {
		cudaSetDevice(dev);
		cudaMemcpy(_devX[dev], xVal, _numBytes, cudaMemcpyHostToDevice);
		cudaMemcpy(_devY[dev], yVal, _numBytes, cudaMemcpyHostToDevice);
	}
}
void corelGPU::run() {
	
	int dev = 0;
    for (int k=0; k<_nentries-1; k++) {
        if(k%1000 == 0) cout << k << endl;
        size_t grid_size = (_nentries-k-1)/_blockSize;
        if((_nentries-k-1) % _blockSize) ++grid_size;

        // launch kernel on both device alternatively
		cudaSetDevice(dev);
		etBim<<<grid_size, _blockSize>>>(_devX[dev], _devY[dev], _devHisto[dev], _nbins, _binWidth, _vMin, _vMax, k, (_nentries-k-1));
		if(dev == 0) {
			dev = 1;
		} else {
			dev = 0;
		}
    }
}
int* corelGPU::getHisto() {
	int* histo = new int[_nbins*_blockSize]; 
	int* tmp = new int[_nbins*_blockSize];
	
	// copy histogram back to host memory space
	cudaSetDevice(0);
	cudaMemcpy(histo, _devHisto[0], _numBytesHisto, cudaMemcpyDeviceToHost);
	cudaSetDevice(1);
	cudaMemcpy(tmp, _devHisto[1], _numBytesHisto, cudaMemcpyDeviceToHost);
	for (int i=0; i<_nbins*_blockSize; i++) {
//		cout << i << " " << histo[i] << " " << tmp[i] << endl;
		histo[i] = histo[i]+tmp[i];
//		cout << histo[i] << endl;
	}
	
	return histo;	
}

size_t corelGPU::getBlockSize() {
	return _blockSize;
}
	
// Main
//
    int main(int argc, char* argv[]) {
    
    double pi = acos(-1.0);

    TApplication* rootapp = new TApplication("Example",&argc,argv);
    TH2* gal = new TH2F("gal", "GPU - Random galaxy phi / theta distribution", 100, 24, 37, 100, 156, 170);

// Open simulated galaxy catalog
//
    TFile *f = new TFile("../Catalogs/Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root");
    TTree *tree = (TTree*)f->Get("bcc");
    bcc *r = new bcc(tree);

    if (r->fChain == 0) return 99;

    Long64_t nentries = r->fChain->GetEntriesFast();

    // define a redshift slice
    float z_min = 0.4;
    float z_max = 0.5;

    Long64_t nbytes=0, nb=0;

    float* ra = new float[nentries];
    float* dec = new float[nentries];

    // Prepare histogram arrays
    // In order to avoid conflict while writing to memory We will have 1 histograms with nbins bins per thread
    //
    int nbins = 50000;
    float ang_min = 0.0;
    float ang_max = 0.2;

    TH1* galgal = new TH1F("galgal", "GPU - Distance galaxy - galaxy", nbins, ang_min, ang_max);

    // Extract ra, dec values from ntuple and copy them into two arrays
    // Select ra and dec only if they belong to the redshift slice
    //
    int nvalues = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
         Long64_t ientry = r->LoadTree(jentry);
        if (ientry < 0) break;

        nb = r->fChain->GetEntry(jentry); 
        float z = r->z;
        if(z>z_min && z<z_max) {
            ra[nvalues] = r->ra*pi/180.0;
            dec[nvalues] = r->dec*pi/180.0;
            nvalues++;
        }
    }

    cout << "Found " << nvalues << " galaxies in z slice" << endl;

	corelGPU* corel = new corelGPU(nvalues, nbins, ang_min, ang_max);
	corel->loadValues(ra, dec);
	corel->run();
	
    int* histo = corel->getHisto();

    // fill root histogram
    for (int i=0; i<nbins; i++) {
        for (int j=0; j<corel->getBlockSize(); j++) {
            int kbin = galgal->GetBinContent(i+1);
            galgal->SetBinContent(i+1,kbin+histo[j*nbins+i]);
        }
    }

    // free memory
    delete ra;
    delete dec;
    delete histo;

	delete corel;

    TFile *h = new TFile("../output/reference_gpu.root","recreate");
    galgal->Write();
    h->Close();

    return 0;
}


