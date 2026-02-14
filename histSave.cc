#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TDirectory.h"


void SaveCorrelationHistograms(
    const char* outFileName,
    const std::vector<TH1F*>& omegaCor,                     //required
    const std::vector<TH1F*>& antiOmegaCor,                  //required
    const std::vector<TH1F*>& backgroundOmegaSum = {},       // optional
    const std::vector<TH1F*>& backgroundAntiOmegaSum = {},   // optional
    const std::vector<TH1F*>& backgroundHighOmega = {},      // optional
    const std::vector<TH1F*>& backgroundHighAntiOmega = {},  // optional
    const std::vector<TH1F*>& backgroundLowOmega = {},       // optional
    const std::vector<TH1F*>& backgroundLowAntiOmega = {}    // optional
) {
    TFile out(outFileName, "RECREATE");


    auto writeVector = [](const std::vector<TH1F*>& vec, const char* prefix) {
        for (size_t i = 0; i < vec.size(); ++i) {
            TH1F* h = vec[i];
            if (!h) continue;
            TH1F* cloneH = (TH1F*)h->Clone(Form("%s_%zu", prefix, i));
            cloneH->SetDirectory(gDirectory);
            cloneH->Write();
            delete cloneH;
        }
    };

    // ========== correlation ==========
    out.mkdir("correlation");
    out.cd("correlation");
    writeVector(omegaCor, "omegaCor");
    writeVector(antiOmegaCor, "antiOmegaCor");

    // ========== background: sum ==========
    if (!backgroundOmegaSum.empty() || !backgroundAntiOmegaSum.empty()) {
        out.cd(); // root
        out.mkdir("background");
        out.cd("background");
        out.mkdir("background/sum");
        out.cd("background/sum");
        writeVector(backgroundOmegaSum, "backgroundOmegaSum");
        writeVector(backgroundAntiOmegaSum, "backgroundAntiOmegaSum");
    }

    // ========== background: high ==========
    if (!backgroundHighOmega.empty() || !backgroundHighAntiOmega.empty()) {
        out.cd(); // root
        // upewnij sie, ze katalog background istnieje
        out.cd();
        out.cd("background");
        out.mkdir("background/high");
        out.cd("background/high");
        writeVector(backgroundHighOmega, "backgroundHighOmega");
        writeVector(backgroundHighAntiOmega, "backgroundHighAntiOmega");
    }

    // ========== background: low ==========
    if (!backgroundLowOmega.empty() || !backgroundLowAntiOmega.empty()) {
        out.cd(); // root
        out.cd("background");
        out.mkdir("background/low");
        out.cd("background/low");
        writeVector(backgroundLowOmega, "backgroundLowOmega");
        writeVector(backgroundLowAntiOmega, "backgroundLowAntiOmega");
    }

    out.Close();
    std::cout << "Zapisano histogramy do pliku: " << outFileName << std::endl;
}

void savekStarDistributions(const char* outFileName,
    const std::vector<TH1F*>& omegakStar,                     //required
    const std::vector<TH1F*>& antiOmegakStar,                  //required
    const std::vector<TH1F*>& backgroundOmegaSum = {},       // optional
    const std::vector<TH1F*>& backgroundAntiOmegaSum = {},   // optional
    const std::vector<TH1F*>& backgroundHighOmega = {},      // optional
    const std::vector<TH1F*>& backgroundHighAntiOmega = {},  // optional
    const std::vector<TH1F*>& backgroundLowOmega = {},       // optional
    const std::vector<TH1F*>& backgroundLowAntiOmega = {}    // optional
) 
{
    TFile out(outFileName, "RECREATE");


    auto writeVector = [](const std::vector<TH1F*>& vec, const char* prefix) {
        for (size_t i = 0; i < vec.size(); ++i) {
            TH1F* h = vec[i];
            if (!h) continue;
            TH1F* cloneH = (TH1F*)h->Clone(Form("%s_%zu", prefix, i));
            cloneH->SetDirectory(gDirectory);
            cloneH->Write();
            delete cloneH;
        }
    };
    std::cout << "Zapisuje histogramy k* do pliku: " << outFileName << std::endl;

    // ========== correlation ==========
    out.mkdir("signal");
    out.cd("signal");
    writeVector(omegakStar, "omegakStar");
    writeVector(antiOmegakStar, "antiOmegakStar");
    std::cout << "Zapisano histogramy k* do pliku: " << outFileName << std::endl;

     // ========== background: sum ==========

    // ========== background: sum ==========
    if (!backgroundOmegaSum.empty() || !backgroundAntiOmegaSum.empty()) {
        out.cd(); // root
        out.mkdir("background");
        out.cd("background");
        out.mkdir("background/sum");
        out.cd("background/sum");
        writeVector(backgroundOmegaSum, "backgroundOmegaSum");
        writeVector(backgroundAntiOmegaSum, "backgroundAntiOmegaSum");
        std::cout << "Zapisano bkg sum k* do pliku: " << outFileName << std::endl;
    }

    // ========== background: high ==========
    if (!backgroundHighOmega.empty() || !backgroundHighAntiOmega.empty()) {
        out.cd(); // root
        // upewnij sie, ze katalog background istnieje
        out.cd();
        out.cd("background");
        out.mkdir("background/high");
        out.cd("background/high");
        writeVector(backgroundHighOmega, "backgroundHighOmega");
        writeVector(backgroundHighAntiOmega, "backgroundHighAntiOmega");
    }

    // ========== background: low ==========
    if (!backgroundLowOmega.empty() || !backgroundLowAntiOmega.empty()) {
        out.cd(); // root
        out.cd("background");
        out.mkdir("background/low");
        out.cd("background/low");
        writeVector(backgroundLowOmega, "backgroundLowOmega");
        writeVector(backgroundLowAntiOmega, "backgroundLowAntiOmega");
    }
    out.Close();
    std::cout << "Zapisano histogramy kstar do pliku: " << outFileName << std::endl;
}
