#define bcc_cxx

#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH2F.h"

#include "bcc.h"


using namespace std;

double angDist(double ra_1, double dec_1, double ra_2, double dec_2) {
//
// Return angular distance between two points in spherical coordinate 
//
    double theta;
    theta = atan(sqrt(cos(dec_2)*cos(dec_2)*sin(ra_2-ra_1)*sin(ra_2-ra_1) + 
                                        (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))*(cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))) /
                                    (sin(dec_1)*sin(dec_2) + cos(dec_1)*cos(dec_2)*cos(ra_2-ra_1)));

    return theta;
}
// Main
//
    int main(int argc, char* argv[]) {
    
    double pi = acos(-1.0);
    int count = 0;

    TApplication* rootapp = new TApplication("Example",&argc,argv);

// Open simulated galaxy catalog
//
    TFile *f = new TFile("/sps/lsst/dev/boutigny/Catalogs/Aardvark/Catalog_v1.0/truth_oscillationcorrected_unrotated/Aardvark_v1.0c_truth.190.root");
    TTree *tree = (TTree*)f->Get("bcc");
    bcc *r = new bcc(tree);

    if (r->fChain == 0) return 99;

    Long64_t nentries = r->fChain->GetEntriesFast();
    Long64_t nbytes=0, nb=0;

    int nvalues = 10000;
    int nbins = 50000;

    double *ra = new double[nvalues];
    double *dec = new double[nvalues];

    for(Long64_t jentry=0; jentry<nvalues; jentry++) {
        Long64_t ientry = r->LoadTree(jentry);
        if (ientry < 0) break;
        nb = r->fChain->GetEntry(jentry);
        ra[jentry] = r->ra*pi/180.0;
        dec[jentry] = r->dec*pi/180.0;
    }

// Compute distances between all pairs of simulated galaxies and store results in galgal histogram
//
    TH1* galgal = new TH1F("galgal", "Distance galaxy - galaxy", nbins, 0.0, 0.2);

    for (int j=0; j<nvalues-1; j++) {
        double ra_1 = ra[j]; double dec_1=dec[j];
        if(j%100 == 0) cout << j << endl;
        for(int k=j+1; k<nvalues; k++) {
            double ra_2 = ra[k]; double dec_2=dec[k];
            double angle = angDist(ra_1, dec_1, ra_2, dec_2);
            galgal->Fill(angle);
        }
    }

    TFile *h = new TFile("reference_no_gpu.root","recreate");
    galgal->Write();
    h->Close();

    return 0;

}


