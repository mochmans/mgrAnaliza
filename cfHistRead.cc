#include <TFile.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <TKey.h>
#include <vector>
#include <string>
#include <iostream>

struct HistogramGroups {
    std::vector<TH1F*> omegaCor;
    std::vector<TH1F*> antiOmegaCor;

    std::vector<TH1F*> backgroundOmegaSum;
    std::vector<TH1F*> backgroundAntiOmegaSum;

    std::vector<TH1F*> backgroundHighOmega;
    std::vector<TH1F*> backgroundHighAntiOmega;

    std::vector<TH1F*> backgroundLowOmega;
    std::vector<TH1F*> backgroundLowAntiOmega;
};

std::vector<TString>getNamesFromDirectory(TFile* file,
                       const char* dirName,
                       const char* prefix)
{
    TDirectory* dir = (TDirectory*)file->Get(dirName);
    std::vector<TString> names;
    if (!dir) {
        std::cout << "Brak katalogu: " << dirName << std::endl;
        TString empty = "";
        std::vector<TString> emptyVec;
        emptyVec.push_back(empty);
        return emptyVec;
    }

    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        std::cout << "Znaleziono obiekt: " << name << std::endl;

        if (!name.BeginsWith(prefix))
            continue;

        names.push_back(name);

    }
    return names;
}

HistogramGroups loadCorrelationHistograms(const char* fileName)

{
    HistogramGroups hists;

    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Nie mozna otworzyc pliku: " << fileName << std::endl;
        return hists;
    }

    // ===== correlation =====
    std::vector<TString> omegaCorNames = getNamesFromDirectory(file, "correlation", "omegaCor");
    std::vector<TString> antiOmegaCorNames = getNamesFromDirectory(file, "correlation", "antiOmegaCor");

    // ===== background sum =====
    std::vector<TString> backgroundOmegaSumNames = getNamesFromDirectory(file, "background/sum", "backgroundOmegaSum");
    std::vector<TString> backgroundAntiOmegaSumNames = getNamesFromDirectory(file, "background/sum", "backgroundAntiOmegaSum");

    // ===== background high =====
    std::vector<TString> backgroundHighOmegaNames = getNamesFromDirectory(file, "background/high", "backgroundHighOmega");
    std::vector<TString> backgroundHighAntiOmegaNames = getNamesFromDirectory(file, "background/high", "backgroundHighAntiOmega");

    // ===== background low =====
    std::vector<TString> backgroundLowOmegaNames = getNamesFromDirectory(file, "background/low", "backgroundLowOmega");
    std::vector<TString> backgroundLowAntiOmegaNames = getNamesFromDirectory(file, "background/low", "backgroundLowAntiOmega");

    TDirectory* corDir = (TDirectory*)file->Get("correlation");
    if(corDir) {    
        if(omegaCorNames[0] != "") {
            for (int i = 0; i < omegaCorNames.size(); i++) {
                TH1F* hist = (TH1F*)corDir->Get(omegaCorNames[i]);
                hists.omegaCor.push_back(hist);
            }
        }
        
        if(antiOmegaCorNames[0] != "") {
            for (int i = 0; i < antiOmegaCorNames.size(); i++) {
                TH1F* hist = (TH1F*)corDir->Get(antiOmegaCorNames[i]);
                hists.antiOmegaCor.push_back(hist);
            }
        }   
    }

    TDirectory* backgroundSumDir = (TDirectory*)file->Get("background/sum");
    if(backgroundSumDir) {
        if(backgroundOmegaSumNames[0] != "") {
            for (int i = 0; i < backgroundOmegaSumNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundSumDir->Get(backgroundOmegaSumNames[i]);
                hists.backgroundOmegaSum.push_back(hist);
            }
        }
        if(backgroundAntiOmegaSumNames[0] != "") {
            for (int i = 0; i < backgroundAntiOmegaSumNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundSumDir->Get(backgroundAntiOmegaSumNames[i]);
                hists.backgroundAntiOmegaSum.push_back(hist);
            }
        }
    }

    TDirectory* backgroundHighDir = (TDirectory*)file->Get("background/high");
    if(backgroundHighDir) {
        if(backgroundHighOmegaNames[0] != "") {
            for (int i = 0; i < backgroundHighOmegaNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundHighDir->Get(backgroundHighOmegaNames[i]);
                hists.backgroundHighOmega.push_back(hist);
            }
        }
        if(backgroundHighAntiOmegaNames[0] != "") {
            for (int i = 0; i < backgroundHighAntiOmegaNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundHighDir->Get(backgroundHighAntiOmegaNames[i]);
                hists.backgroundHighAntiOmega.push_back(hist);
            }
        }
    }

    TDirectory* backgroundLowDir = (TDirectory*)file->Get("background/low");
    if(backgroundLowDir) {
        if(backgroundLowOmegaNames[0] != "") {
            for (int i = 0; i < backgroundLowOmegaNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundLowDir->Get(backgroundLowOmegaNames[i]);
                hists.backgroundLowOmega.push_back(hist);
            }
        }
        if(backgroundLowAntiOmegaNames[0] != "") {
            for (int i = 0; i < backgroundLowAntiOmegaNames.size(); i++) {
                TH1F* hist = (TH1F*)backgroundLowDir->Get(backgroundLowAntiOmegaNames[i]);
                hists.backgroundLowAntiOmega.push_back(hist);
            }
        }
    }

    return hists;

}