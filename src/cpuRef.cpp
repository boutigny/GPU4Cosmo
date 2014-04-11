#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#define bcc_cxx
#include "bcc.h"


using namespace std;

//
// Return angular distance between two points in spherical coordinates.
// ra_ and dec_ are right ascension and declination in radians for the two points.
//
float AngDistance(float ra_1, float dec_1, float ra_2, float dec_2)
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
    cout << "Usage: " << progName << " <inputFile> <outpuFile> <number of bins>" << endl;
}


//
// Main
//
int main(int argc, const char* argv[])
{
    // Parse command line
    if (argc < 4) {
        Usage(argv[0]);
        return 1;
    }
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];
    const int kNumBins = atoi(argv[3]); 
    

    // Open simulated galaxy catalog
    // const char* inputFileName = "../Catalogs/Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root";
    TFile* inFile = TFile::Open(inputFileName);
    if (inFile == 0) {
        cout << "could not open file '" << inputFileName << "'" << endl;
        return 2;
    }

    // Load the spherical coordinates of the galaxies from the input ROOT
    // file
    TTree* tree = (TTree*)inFile->Get("bcc");
    bcc* r = new bcc(tree);
    if (r->fChain == 0) {
        return 99;
    }

    int kMaxGalaxies = 10000;
    float* ra = new float[kMaxGalaxies];
    float* dec = new float[kMaxGalaxies];
    Long64_t numGalaxies;
    for (numGalaxies=0; numGalaxies < kMaxGalaxies; numGalaxies++) {
        if (r->LoadTree(numGalaxies) < 0) {
            break;
        }
        if (r->GetEntry(numGalaxies) == 0) {
            break;
        }
        ra[numGalaxies] = DegToRad(r->ra);
        dec[numGalaxies] = DegToRad(r->dec);
    }

    if (numGalaxies < 2) {
        cout << "Found " << numGalaxies << " galaxies in file "
             << inputFileName << endl
             << "Aborting execution" << endl;
        return 99;
    }


    // Compute distances between all pairs of simulated galaxies and store
    // results in 'galgal' histogram in ROOT
    cout << "Processing " << numGalaxies << " galaxies" << endl;

    const double kBinLow = 0.0;
    const double kBinUp = 0.2;
    TH1* histo = new TH1F("galgal", "Distance galaxy - galaxy", kNumBins, kBinLow,
                          kBinUp);
    for (int j=0; j < numGalaxies-1; j++) {
        for (int k=j+1; k < numGalaxies; k++) {
            histo->Fill(AngDistance(ra[j], dec[j], ra[k], dec[k]));
        }

        if (j%100 == 0) {
            cout << j << endl;
        }
    }

    // Close files
    // Prepare output file
    // TFile* outFile = new TFile("../output/reference_no_gpu.root", "recreate");
    TFile* outFile = TFile::Open(outputFileName, "recreate");
    if (outFile == 0) {
        cout << "could not create file '" << outputFileName << "'" << endl;
        return 2;
    }
    histo->Write();
    outFile->Close();
    inFile->Close();

    return 0;
}
