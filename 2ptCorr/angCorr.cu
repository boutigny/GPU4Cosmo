#define bcc_cxx

#include <iostream>
#include <assert.h>
#include "math.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH2F.h"

#include "healpix_base.h"
#include "healpix_map.h"
#include "planck_rng.h"
#include "pointing.h"
#include "error_handling.h"

#include "bcc.h"

using namespace std;

__device__ float angDist(float ra_1, float dec_1, float ra_2, float dec_2) {
//
// Return angular distance between two points in spherical coordinate 
//
//    float numerator, denominator, theta;
	
//	numerator = sqrt(cos(dec_2)*cos(dec_2)*sin(ra_2-ra_1)*sin(ra_2-ra_1) + 
//                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))*
//                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1)));
//	denominator = (sin(dec_1)*sin(dec_2) + cos(dec_1)*cos(dec_2)*cos(ra_2-ra_1));
	
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

	
//	theta = atan2(numerator,denominator);
	
//    theta = atan(sqrt(cos(dec_2)*cos(dec_2)*sin(ra_2-ra_1)*sin(ra_2-ra_1) + 
//                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))*
//                                (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))) /
//                                (sin(dec_1)*sin(dec_2) + cos(dec_1)*cos(dec_2)*cos(ra_2-ra_1)));

//    return theta;
}

__global__ void etBim(const float *ra, const float *dec, int *histo, const int nbins, const float bin_width, const float ang_min, 
                    const float ang_max, const int index, const int n) {

// compute the angular separation between every object pair identified by their position in the ra and dec arrays//

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
__global__ void etBimBam(const float *ra, const float *dec, const float *ra2, const float *dec2, 
                    int *histo, const int nbins, const float bin_width, const float ang_min, 
                    const float ang_max, const int index, const int n) {

// compute angular separation between all objects in a first set (ra, dec) and their counterpart in a second set (ra2, dec2)
//

    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;

    // There are 2 possibilities to compute the correlation between 2 different field,
    // not completly sure which one is correct.

    if(i < n) {
        float angle = angDist(ra[index], dec[index], ra2[i], dec2[i]);
        if (angle < ang_max) {
            int bin = (int)((angle-ang_min)/bin_width);
            int ibin = threadIdx.x*nbins+bin;
            histo[ibin]++;
        }
    }

/*	if(i < n) {
		float angle = angDist(ra[index], dec[index], ra2[i+index+1], dec2[i+index+1]);
        if (angle < ang_max) {
            int bin = (int)((angle-ang_min)/bin_width);
            int ibin = threadIdx.x*nbins+bin;
            histo[ibin]++;
        }
    }
*/
}

class rndField {

// Generate a random field of galaxy positions (ra, dec) within a spherical tile
// defined acoording to the HEALPix scheme
//
    public:
        rndField(int);
        double *getRaDec();
        void genField(int, int);
        double* getMap();
    private:
        Healpix_Base _base;
        int _nside;
        double* _map;
};

rndField::rndField(int nside) : _nside(nside) {
    
    _base = Healpix_Base(_nside,RING);
}

void rndField::genField(int cell, int values) {

    _map = new double[2*values];

    pointing p = _base.pix2ang(cell);
    double pi = acos(-1.0);
    double th_c = p.theta;
    double ph_c = p.phi;

    double cth_min = cos(th_c-7.0*pi/180.0);
    double cth_max = cos(th_c+7.0*pi/180.0);
    if (cth_min > cth_max) {
        int tmp = cth_max;
        cth_max = cth_min;
        cth_min = tmp;
    }
    double ph_min = ph_c-7.0*pi/180.0;
    double ph_max = ph_c+7.0*pi/180.00;

    //	Use the Planck random number generator available in HEALPix
    planck_rng rng;

    int count = 0;
    while (count < values) {
       double phi = (ph_max-ph_min)*rng.rand_uni()+ph_min;
       double theta = acos( (cth_max-cth_min)*rng.rand_uni()+cth_min );

        pointing pt(theta,phi);
        int pix = _base.ang2pix(pt);
        if (pix != cell) {
            continue;
        } 
        _map[count] = phi;
        _map[values+count+1] = pi/2.0 - theta;
        count++;
            
//        cout << count << " " << theta*180.0/pi << " - " << phi*180.0/pi << " "  << pix <<  endl;
    }
}

double* rndField::getMap() {
    return _map;
}

class corelGPU {
// Compute 2 point correlation function on GPU
//
	public:
		corelGPU(int, int, float, float);
		~corelGPU();
		void loadValues(float*, float*);
		void loadValues(float*, float*, float*, float*);
		void run(int);
		int* getHisto();
		void resetHisto();
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
		
		// _devX and _devY are 2 floating point arrays on device memory containing the values of 
		// the variables to be corelated
		float* _devX;
		float* _devY;
		
		// _devX2 and _devY2 : second set of quantities to be correlated when we
		// we correlate 2 diffferent fields
		float* _devX2;
		float* _devY2;
		
		// _devHisto is an histogram containing the result of the correlation computation (distribution of the correlation values)
		int* _devHisto;
};

corelGPU::corelGPU(int nentries, int nbins, float vMin, float vMax) : _nentries(nentries), _nbins(nbins),
			_vMin(vMin), _vMax(vMax) 
{
	_numBytes = nentries*sizeof(float);
	_numBytesHisto = nbins*_blockSize*sizeof(int);
	
	_devX = 0;
	_devY = 0;
	_devHisto =0;
	
	_devX2 = 0;
	_devY2 = 0;
	
	_binWidth = (vMax-vMin)/(float)nbins;
    
	// allocate device memory
    //
    cudaMalloc((void**)&_devX, _numBytes);
    cudaMalloc((void**)&_devY, _numBytes);

    cudaMalloc((void**)&_devHisto, _numBytesHisto);
    cudaMemset(_devHisto,0,_numBytesHisto);
	
	// check allocation
    //
    assert(_devX != 0 && _devY != 0 && _devHisto != 0);	
}
corelGPU::~corelGPU() {

    cudaFree(_devX);
    cudaFree(_devY);
    cudaFree(_devHisto);
}
void corelGPU::loadValues(float* xVal, float* yVal) {

    // Copy xVal and yVal arrays to device memory
    //
    cudaMemcpy(_devX, xVal, _numBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(_devY, yVal, _numBytes, cudaMemcpyHostToDevice);
}
void corelGPU::loadValues(float* xVal, float* yVal, float* x2Val, float* y2Val) {

	// Allocate device memory for second set of values to be correlated
	//
	if(_devX2 == 0) cudaMalloc((void**)&_devX2, _numBytes);
	if(_devY2 == 0) cudaMalloc((void**)&_devY2, _numBytes);
	assert(_devX2 != 0 &&_devY2 !=0 );
	
    // Copy arrays to device memory
    //
    cudaMemcpy(_devX, xVal, _numBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(_devY, yVal, _numBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(_devX2, x2Val, _numBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(_devY2, y2Val, _numBytes, cudaMemcpyHostToDevice);
}
void corelGPU::run(int nsets) {

	assert(nsets == 1 || nsets == 2);
	
	if(nsets == 1) {
		for (int k=0; k<_nentries-1; k++) {
			if(k%1000 == 0) cout << k << endl;
			size_t grid_size = (_nentries-k-1)/_blockSize;
			if((_nentries-k-1) % _blockSize) ++grid_size;

			// launch kernel
			etBim<<<grid_size, _blockSize>>>(_devX, _devY, _devHisto, _nbins, _binWidth, _vMin, _vMax, k, (_nentries-k-1));
		}
	}
	if(nsets == 2) {
		size_t grid_size = _nentries/_blockSize;
		if(_nentries % _blockSize) ++grid_size;
		for (int k=0; k<_nentries; k++) {
			if(k%1000 ==0) cout << k << endl;
			// launch kernel
			etBimBam<<<grid_size, _blockSize>>>(_devX, _devY, _devX2, _devY2, _devHisto, _nbins, _binWidth, _vMin, _vMax, k, _nentries);
		}

/*		for (int k=0; k<_nentries-1; k++) {
			if(k%1000 == 0) cout << k << endl;
			size_t grid_size = (_nentries-k-1)/_blockSize;
			if((_nentries-k-1) % _blockSize) ++grid_size;

			// launch kernel
			etBimBam<<<grid_size, _blockSize>>>(_devX, _devY, _devX2, _devY2, _devHisto, _nbins, _binWidth, _vMin, _vMax, k, (_nentries-k-1));
		}
*/
	}
		
}
int* corelGPU::getHisto() {
	int* histo = new int[_nbins*_blockSize]; 
	
	// copy histogram back to host memory space
    cudaMemcpy(histo, _devHisto, _numBytesHisto, cudaMemcpyDeviceToHost);
	
	return histo;	
}

void corelGPU::resetHisto() {

    //Reinitialize histogram
    cudaMemset(_devHisto,0,_numBytesHisto);
}

size_t corelGPU::getBlockSize() {
	return _blockSize;
}
	
// Main
//
    int main(int argc, char* argv[]) {
    
    double pi = acos(-1.0);

    TApplication* rootapp = new TApplication("Example",&argc,argv);

// Open simulated galaxy catalog in root format
//
    TFile *f = new TFile("../Catalogs/Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root");
    TTree *tree = (TTree*)f->Get("bcc");
    bcc *r = new bcc(tree);

    if (r->fChain == 0) return 99;

    Long64_t nentries = r->fChain->GetEntriesFast();

    // define a redshift slice
    float z_min = 0.0;
    float z_max = 0.2;

    Long64_t nbytes=0, nb=0;

    float* ra = new float[nentries];
    float* dec = new float[nentries];

    // Prepare histogram arrays
    // In order to avoid conflicts between threads while writing to memory 
	// we use 1 histogram with nbins bins per thread
    //
    int nbins = 50000;
    float ang_min = 0.0;
    float ang_max = 0.2;

    TH1* galgal = new TH1F("galgal", "GPU - Distance galaxy - galaxy", nbins, ang_min, ang_max);
    TH1* zdis = new TH1F("zdis", "z distribution", 200, 0., 2.);

    // Extract ra, dec values from ntuple and copy them into two arrays
    // Select ra and dec only if they belong to the redshift slice
    //
    int nvalues = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
         Long64_t ientry = r->LoadTree(jentry);
        if (ientry < 0) break;

        nb = r->fChain->GetEntry(jentry); 
        float z = r->z;
		zdis->Fill(z);
        if(z>z_min && z<z_max) {
            ra[nvalues] = r->ra*pi/180.0;
            dec[nvalues] = r->dec*pi/180.0;
            nvalues++;
        }
    }

    cout << "Found " << nvalues << " galaxies in z slice" << endl;

    corelGPU* corel = new corelGPU(nvalues, nbins, ang_min, ang_max);
    corel->loadValues(ra, dec);
    corel->run(1);
	
    int* histo = corel->getHisto();

    // fill root histogram
    for (int i=0; i<nbins; i++) {
        for (int j=0; j<corel->getBlockSize(); j++) {
            int kbin = galgal->GetBinContent(i+1);
            galgal->SetBinContent(i+1,kbin+histo[j*nbins+i]);
        }
    }
	
    corel->resetHisto();
	
    // Generate random field of galaxy positions with the same number of
    // galaxies as in the original simulated sample (nvalues)
    //
    int nsides = 3;   // Define tiling granularity in HEALPix
    int cell = 190;   // Tile number (should match the catalog name)
    rndField* field = new rndField(nsides);
    field -> genField(cell, nvalues);
	
    double *map = field->getMap();
	
    // Create histogram to check that theta / phi distribution is as expected
    //
    TH2* gal = new TH2F("gal", "Random galaxy phi / theta distribution", 100, 24, 37, 100, 156, 170);
	
    // compute 2 point correlation function from galaxy positions in the random field
    //
    float* rnd_ra = new float[nvalues];
    float* rnd_dec = new float[nvalues];
	
    for (int i=0; i<nvalues; i++) {
        rnd_ra[i] = map[i]; 
        rnd_dec[i] = map[nvalues+i+1];
        gal->Fill(rnd_dec[i]*180.0/pi,rnd_ra[i]*180.0/pi);
    } 
	
    corel->loadValues(rnd_ra, rnd_dec);
    corel->run(1);
    histo = corel->getHisto();
	
    // fill root histogram
    //
    TH1* rndrnd = new TH1F("rndrnd", "GPU - Distance random - random", nbins, ang_min, ang_max);

    for (int i=0; i<nbins; i++) {
        for (int j=0; j<corel->getBlockSize(); j++) {
            int kbin = rndrnd->GetBinContent(i+1);
            rndrnd->SetBinContent(i+1,kbin+histo[j*nbins+i]);
        }
    }
	
    corel->resetHisto();
	
    // compute 2 point correlation function between galaxy and random field
    //
	
    corel->loadValues(ra, dec, rnd_ra, rnd_dec);
    corel->run(2);
    histo = corel->getHisto();

    // fill root histogram
    //
    TH1* galrnd = new TH1F("galrnd", "GPU - Distance galaxy - random", nbins, ang_min, ang_max);

    for (int i=0; i<nbins; i++) {
        for (int j=0; j<corel->getBlockSize(); j++) {
            int kbin = galrnd->GetBinContent(i+1);
            galrnd->SetBinContent(i+1,kbin+histo[j*nbins+i]);
        }
    }
	
    // Compute random field corrected correlation function
    //
    TH1* corr = new TH1F("corr", "Corrected angular correlation function", nbins, ang_min, ang_max);
    //Normalization
    //
    float ggnorm = (nvalues*nvalues - nvalues)/2.0;
    float rrnorm = ggnorm;
    float grnorm = nvalues*nvalues;
    corr -> Add(galgal, galrnd, 1.0/ggnorm, -2.0/grnorm);
    corr -> Add(corr, rndrnd, 1.0, 1.0/rrnorm);
    corr -> Divide(corr, rndrnd, 1.0, 1.0/rrnorm);

    // free memory
    delete ra;
    delete dec;
    delete rnd_ra;
    delete rnd_dec;
    delete histo;

    delete corel;

    TFile *h = new TFile("../output/gpu.root","recreate");
    galgal->Write();
    gal->Write();
    zdis->Write();
    rndrnd->Write();
    galrnd->Write();
    corr->Write();
    h->Close();

    return 0;
}


