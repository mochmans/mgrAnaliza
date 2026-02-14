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
#include "backgroundCorrelation.cc"
#include "histSave.cc"
#include "histRead.cc"

HistogramGroups kStarDist(TString rootFileName, bool saveHists = false){

    // rootFileName = "mediumDataOmegaARfull.root";
    // BUGGED ZAPIS ????????? ////
    TFile* f = new TFile(rootFileName,"open");

    HistogramGroups hists;

    hists.omegaCor = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended");
    hists.antiOmegaCor = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended_antyomega");

    // hists.omegaCor[0]->Draw();
    // for(size_t i = 1; i < hists.omegaCor.size(); i++) {
    //     std::cout << "Adding histogram: " << hists.omegaCor[i]->GetName() << std::endl;
    //     hists.omegaCor[i]->Draw("SAME");
    // }
    

    hists.backgroundOmegaSum = kstarDistributionsSum(f,"femto-universe-pair-task-track-cascade-extended_omega_bkg_low","femto-universe-pair-task-track-cascade-extended_omega_bkg_high");
    hists.backgroundAntiOmegaSum = kstarDistributionsSum(f,"femto-universe-pair-task-track-cascade-extended_antyomega_bkg_low","femto-universe-pair-task-track-cascade-extended_antyomega_bkg_high");

    hists.backgroundLowOmega = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended_omega_bkg_low");
    hists.backgroundLowAntiOmega = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended_antyomega_bkg_low");

    hists.backgroundHighOmega = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended_omega_bkg_high");   
    hists.backgroundHighAntiOmega = kstarDistributions(f,"femto-universe-pair-task-track-cascade-extended_antyomega_bkg_high");

    

    std::cout << "Loaded k* distributions from file: " << rootFileName << std::endl;

    if(saveHists) {
        savekStarDistributions(
            rootFileName.ReplaceAll(".root","_kStarDists.root"),
            hists.omegaCor,
            hists.antiOmegaCor,
            hists.backgroundOmegaSum,
            hists.backgroundAntiOmegaSum,
            hists.backgroundHighOmega,
            hists.backgroundHighAntiOmega,
            hists.backgroundLowOmega,
            hists.backgroundLowAntiOmega
        );
    }
    f->Close();
    return hists;
}

void kStarMerge(std::vector<TString> rootFileNames, const char* outFileName){
    std::vector<HistogramGroups> allHists;
    
    for(const auto& fileName : rootFileNames) {
        allHists.push_back(kStarDist(fileName, false));
    }

    std::vector<TH1F*> mergedOmegakStar;
    std::vector<TH1F*> mergedAntiOmegakStar;
    std::vector<TH1F*> mergedBackgroundOmegaSum;
    std::vector<TH1F*> mergedBackgroundAntiOmegaSum;
    std::vector<TH1F*> mergedBackgroundHighOmega;
    std::vector<TH1F*> mergedBackgroundHighAntiOmega;
    std::vector<TH1F*> mergedBackgroundLowOmega;
    std::vector<TH1F*> mergedBackgroundLowAntiOmega;

    TCanvas* c = new TCanvas("Merged k* distributions", "Merged k* distributions", 800, 600);
    c->Divide(2,1);
    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> omegaCorToMerge;
        for(const auto& hists : allHists) {
            std::cout << "Adding histogram: " << hists.omegaCor[i]->GetName() << std::endl;
            omegaCorToMerge.push_back((TH1F*)hists.omegaCor[i]->Clone(Form("temp_%d", i)));
            
        }
        TH1F* mergedOmegaCorHisto;
        mergedOmegaCorHisto = (TH1F*)omegaCorToMerge[0]->Clone("mergedOmegakStar");
        for(int j = 1; j < omegaCorToMerge.size(); j++) {
            mergedOmegaCorHisto->Add(omegaCorToMerge[j]);
        }

        mergedOmegakStar.push_back(mergedOmegaCorHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> antiOmegaCorToMerge;
        for(const auto& hists : allHists) {
            antiOmegaCorToMerge.push_back(hists.antiOmegaCor[i]);
        }
        TH1F* mergedAntiOmegaCorHisto = (TH1F*)antiOmegaCorToMerge[0]->Clone(Form("mergedAntiOmegakStar_%d", i));
        for(int j = 1; j < antiOmegaCorToMerge.size(); j++) {
            mergedAntiOmegaCorHisto->Add(antiOmegaCorToMerge[j]);
        }
        mergedAntiOmegakStar.push_back(mergedAntiOmegaCorHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundOmegaSumToMerge;
        for(const auto& hists : allHists) {
            backgroundOmegaSumToMerge.push_back(hists.backgroundOmegaSum[i]);
        }
        TH1F* mergedBackgroundOmegaSumHisto = (TH1F*)backgroundOmegaSumToMerge[0]->Clone(Form("mergedBackgroundOmegaSum_%d", i));
        for(int j = 1; j < backgroundOmegaSumToMerge.size(); j++) {
            mergedBackgroundOmegaSumHisto->Add(backgroundOmegaSumToMerge[j]);
        }
        mergedBackgroundOmegaSum.push_back(mergedBackgroundOmegaSumHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundAntiOmegaSumToMerge;
        for(const auto& hists : allHists) {
            backgroundAntiOmegaSumToMerge.push_back(hists.backgroundAntiOmegaSum[i]);   
        }
        TH1F* mergedBackgroundAntiOmegaSumHisto = (TH1F*)backgroundAntiOmegaSumToMerge[0]->Clone(Form("mergedBackgroundAntiOmegaSum_%d", i));
        for(int j = 1; j < backgroundAntiOmegaSumToMerge.size(); j++) {
            mergedBackgroundAntiOmegaSumHisto->Add(backgroundAntiOmegaSumToMerge[j]);
        }
        mergedBackgroundAntiOmegaSum.push_back(mergedBackgroundAntiOmegaSumHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundHighOmegaToMerge;
        for(const auto& hists : allHists) {
            backgroundHighOmegaToMerge.push_back(hists.backgroundHighOmega[i]);
        }
        TH1F* mergedBackgroundHighOmegaHisto = (TH1F*)backgroundHighOmegaToMerge[0]->Clone(Form("mergedBackgroundHighOmega_%d", i));
        for(int j = 1; j < backgroundHighOmegaToMerge.size(); j++) {
            mergedBackgroundHighOmegaHisto->Add(backgroundHighOmegaToMerge[j]);
        }
        mergedBackgroundHighOmega.push_back(mergedBackgroundHighOmegaHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundHighAntiOmegaToMerge;
        for(const auto& hists : allHists) {
            backgroundHighAntiOmegaToMerge.push_back(hists.backgroundHighAntiOmega[i]);
        }
        TH1F* mergedBackgroundHighAntiOmegaHisto = (TH1F*)backgroundHighAntiOmegaToMerge[0]->Clone(Form("mergedBackgroundHighAntiOmega_%d", i));
        for(int j = 1; j < backgroundHighAntiOmegaToMerge.size(); j++) {
            mergedBackgroundHighAntiOmegaHisto->Add(backgroundHighAntiOmegaToMerge[j]);
        }
        mergedBackgroundHighAntiOmega.push_back(mergedBackgroundHighAntiOmegaHisto);
    }

    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundLowOmegaToMerge;
        for(const auto& hists : allHists) {
            backgroundLowOmegaToMerge.push_back(hists.backgroundLowOmega[i]);
        }
        TH1F* mergedBackgroundLowOmegaHisto = (TH1F*)backgroundLowOmegaToMerge[0]->Clone(Form("mergedBackgroundLowOmega_%d", i));
        for(int j = 1; j < backgroundLowOmegaToMerge.size(); j++) {
            mergedBackgroundLowOmegaHisto->Add(backgroundLowOmegaToMerge[j]);
        }
        mergedBackgroundLowOmega.push_back(mergedBackgroundLowOmegaHisto);
    }
    
    for(int i = 0; i < 12; i++) {
        std::vector<TH1F*> backgroundLowAntiOmegaToMerge;
        for(const auto& hists : allHists) {
            backgroundLowAntiOmegaToMerge.push_back(hists.backgroundLowAntiOmega[i]);
        }
        TH1F* mergedBackgroundLowAntiOmegaHisto = (TH1F*)backgroundLowAntiOmegaToMerge[0]->Clone(Form("mergedBackgroundLowAntiOmega_%d", i));
        for(int j = 1; j < backgroundLowAntiOmegaToMerge.size(); j++) {
            mergedBackgroundLowAntiOmegaHisto->Add(backgroundLowAntiOmegaToMerge[j]);
        }
        mergedBackgroundLowAntiOmega.push_back(mergedBackgroundLowAntiOmegaHisto);
    }

    savekStarDistributions(
            outFileName,
            mergedOmegakStar,
            mergedAntiOmegakStar,
            mergedBackgroundOmegaSum,
            mergedBackgroundAntiOmegaSum,
            mergedBackgroundHighOmega,
            mergedBackgroundHighAntiOmega,
            mergedBackgroundLowOmega,
            mergedBackgroundLowAntiOmega
        );

}

void mainTest(){
    std::vector<TString> rootFileNames = {
        "mediumDataOmegaARfull.root"
        ,"zzh_medium.root"
        ,"zzo_medium.root"
        // ,"zzm_medium.root"
    };

    kStarMerge(rootFileNames, "kStarMergedNozzm.root");
}