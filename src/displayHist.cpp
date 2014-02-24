/*
 * ROOT macro to display histograms
 *
 */

void displayGalGalHisto(const char* fileName="reference_no_gpu.root")
{
    displayHisto(fileName, "galgal");
}


void displayHisto(const char* fileName, const char* histoName)
{
    TFile* f = TFile::Open(fileName);
    if (f == 0) {
        return;
    }

    TH1F* histo = (TH1F*)f->Get(histoName);
    histo->Draw();
}

