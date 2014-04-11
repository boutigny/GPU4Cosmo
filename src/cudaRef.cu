#include <iostream>
#include <assert.h>
#include "math.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH2F.h"

#define bcc_cxx
#include "bcc.h"

using namespace std;

//
// Return angular distance between two points in spherical coordinate
//
__device__ float angDist(float ra_1, float dec_1, float ra_2, float dec_2)
{
    float sindec1 = sinf(dec_1);
    float cosdec1 = cosf(dec_1);
    float sindec2 = sinf(dec_2);
    float cosdec2 = cosf(dec_2);
    float cosra2_ra1 = cosf(ra_2-ra_1);
    float sinra2_ra1 = sinf(ra_2-ra_1);

    float aux = (cosdec1 * sindec2) - (sindec1 * cosdec2 * cosra2_ra1);
    float num = (cosdec2 * cosdec2 * sinra2_ra1 * sinra2_ra1) + (aux * aux);
    float den = (sindec1 * sindec2) + (cosdec1 * cosdec2 * cosra2_ra1);

    return atan2f(sqrtf(num), den);
}


__global__ void etBim(const float* ra, const float* dec,
                      int* histo, int nbins, float bin_width,
                      float ang_min,float ang_max,
                      int index, int n)
{
    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < n) {
        float angle = angDist(ra[index], dec[index], ra[i+index+1], dec[i+index+1]);
        if (angle < ang_max) {
            int bin = (int)((angle - ang_min)/bin_width);
            int ibin = threadIdx.x*nbins + bin;
            histo[ibin]++;
        }
    }
}


class corelGPU {
public:
    corelGPU(int, int, float, float, int, int);
    ~corelGPU();
    void loadValues(float*, float*);
    void run();
    int* getHisto();
    int getBlockSize();

private:
    int _nentries;
    int _nbins;
    float _vMin;
    float _vMax;
    static const int _blockSize = 1024;
    size_t _numBytes;
    size_t _numBytesHisto;
    float _binWidth;
    int _dev_1, _dev_2;

    // _devX and _devY are 2 floating point arrays on device memory containing
    // the values of the variables to be corelated.
    // There will be 1 array per device - The considered device are specified by _dev_1 and _dev_2
    float** _devX;
    float** _devY;

    // _devHisto is an histogram containing the result of the correlation computation (distribution of the correlation values)
    // 1 array per device
    int** _devHisto;
};


corelGPU::corelGPU(int nentries, int nbins, float vMin,
                   float vMax, int dev_1, int dev_2) : _nentries(nentries), _nbins(nbins),
    _vMin(vMin), _vMax(vMax), _dev_1(dev_1), _dev_2(dev_2)
{
    _numBytes = nentries*sizeof(float);
    _numBytesHisto = nbins*_blockSize*sizeof(int);

    _devX = new float*[_dev_2-_dev_1+1];
    _devY = new float*[_dev_2-_dev_1+1];
    _devHisto = new int*[_dev_2-_dev_1+1];

    _binWidth = (vMax-vMin)/(float)nbins;

    // allocate device memory on each GPU
    //
    for (int dev=_dev_1; dev<_dev_2+1; dev++) {
        _devX[dev] = 0;
        _devY[dev] = 0;
        _devHisto[dev] = 0;
        cudaSetDevice(dev);
        cudaMalloc((void**)&_devX[dev], _numBytes);
        cudaMalloc((void**)&_devY[dev], _numBytes);
        cudaMalloc((void**)&_devHisto[dev], _numBytesHisto);
        cudaMemset(_devHisto[dev], 0, _numBytesHisto);

        // check allocation
        //
        assert(_devX[dev] != 0 && _devY[dev] != 0 && _devHisto[dev] != 0);
    }
}


corelGPU::~corelGPU()
{

    for (int dev=_dev_1; dev<_dev_2+1; dev++) {
        cudaSetDevice(dev);
        cudaFree(_devX[dev]);
        cudaFree(_devY[dev]);
        cudaFree(_devHisto[dev]);
    }
}


void corelGPU::loadValues(float* xVal, float* yVal)
{

    // Copy xVal and yVal arrays to device memory
    //
    for (int dev=_dev_1; dev < _dev_2+1; dev++) {
        cudaSetDevice(dev);
        cudaMemcpy(_devX[dev], xVal, _numBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(_devY[dev], yVal, _numBytes, cudaMemcpyHostToDevice);
    }
}


void corelGPU::run()
{

    int dev = _dev_1;
    for (int k=0; k<_nentries-1; k++) {
        if (k%1000 == 0) {
            cout << k << endl;
        }
        size_t grid_size = (_nentries-k-1)/_blockSize;
        if ((_nentries-k-1) % _blockSize) {
            ++grid_size;
        }

        // launch kernel on both device alternatively
        cudaSetDevice(dev);
        etBim<<<grid_size, _blockSize>>>(_devX[dev], _devY[dev], _devHisto[dev],
                                         _nbins, _binWidth, _vMin, _vMax, k,
                                         (_nentries-k-1));
        dev++;
        if (dev > _dev_2) {
            dev = _dev_1;
        }
    }
}


int* corelGPU::getHisto()
{
    int* histo = new int[_nbins*_blockSize];
    int* tmp = new int[_nbins*_blockSize];

    // copy histogram back to host memory space
    for (int dev=_dev_1; dev<_dev_2+1; dev++) {
        cudaSetDevice(dev);
        cudaMemcpy(tmp, _devHisto[dev], _numBytesHisto, cudaMemcpyDeviceToHost);
        for (int i=0; i<_nbins*_blockSize; i++) {
            //      cout << i << " " << histo[i] << " " << tmp[i] << endl;
            if (dev == _dev_1) {
                histo[i] = tmp[i];
            }
            else {
                histo[i] = histo[i]+tmp[i];
            }
        }
    }

    return histo;
}


int corelGPU::getBlockSize()
{
    return _blockSize;
}


//
// Converts a value in degrees to radians
//
inline float DegToRad(float deg)
{
    const float PI = acosf(-1.0);
    return deg * (PI/180.0f);
}


//
// Show this program usage
//
void Usage(const char* path)
{
    const char* slash = strrchr(path, '/');
    const char* progName = (slash != NULL) ? ++slash : path;
    cout << "Usage: " << progName << " <inputFile> <outpuFile> <number_of_bins_in_histogram> <first_device> <last_device> (start at 0) <number of galaxies>" << endl;
}


// Main
//
int main(int argc, const char* argv[])
{
    // Parse command line
    if (argc < 7) {
        Usage(argv[0]);
        return 1;
    }
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];
    
    int nbins = atoi(argv[3]);
    int dev_1 = atoi(argv[4]);
    int dev_2 = atoi(argv[5]);
    int kMaxGalaxies = atoi(argv[6]);
    
    cout << "Using device " << dev_1 << " to " << dev_2 << endl;

    // Open simulated galaxy catalog
    //
    // TFile *f = TFile::Open("../Catalogs/Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root");
    TFile* inFile = TFile::Open(inputFileName);
    TTree* tree = (TTree*)inFile->Get("bcc");
    bcc* r = new bcc(tree);
    if (r->fChain == 0) {
        return 99;
    }

    // Define a redshift slice
    float z_min = 0.4;
    float z_max = 0.5;

    // Prepare histogram arrays
    // In order to avoid conflict while writing to memory We will have 1
    // histogram with nbins bins per thread
    float ang_min = 0.0;
    float ang_max = 0.2;

    // Extract ra, dec values from ntuple and copy them into two arrays
    // Select ra and dec only if they belong to the redshift slice
    //
    Long64_t nentries = r->fChain->GetEntriesFast();
    float* ra = new float[nentries];
    float* dec = new float[nentries];
    int nvalues = 0;
    for (Long64_t jentry=0; jentry < nentries; jentry++) {
        if (r->LoadTree(jentry) < 0) {
            break;
        }

        if (r->GetEntry(jentry) == 0) {
            break;
        }

        float z = r->z;
        if (true || (z > z_min && z < z_max)) {
            ra[nvalues] = DegToRad(r->ra);
            dec[nvalues] = DegToRad(r->dec);
            nvalues++;
        }

        if (nvalues == kMaxGalaxies) break;
    }

    // Do the correlation computation
    cout << "Found " << nvalues << " galaxies in z slice" << endl;
    corelGPU* corel = new corelGPU(nvalues, nbins, ang_min, ang_max, dev_1, dev_2);
    corel->loadValues(ra, dec);
    corel->run();

    // Fill ROOT histogram
    TH1* galgal = new TH1F("galgal", "GPU - Distance galaxy - galaxy", nbins,
                           ang_min, ang_max);
    int* histo = corel->getHisto();
    for (int i=0; i < nbins; i++) {
        for (int j=0; j < corel->getBlockSize(); j++) {
            int kbin = galgal->GetBinContent(i + 1);
            galgal->SetBinContent(i+1, kbin+histo[j*nbins + i]);
        }
    }
    galgal->ResetStats();

    // free memory
    delete ra;
    delete dec;
    delete histo;
    delete corel;

    // Write ouput histogram to ROOT file
    // TFile *h = TFile::Open("../output/reference_gpu.root", "recreate");
    TFile* outFile = TFile::Open(outputFileName, "recreate");
    if (outFile == 0) {
        cout << "could not create file '" << outputFileName << "'" << endl;
        return 2;
    }
    galgal->Write();
    outFile->Close();
    inFile->Close();

    return 0;
}


