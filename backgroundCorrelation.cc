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

TH1F *createCFHistogram(TH1F* hNumerator, TH1F* hDenominator, const char* name, Int_t rebinFactor=1){
    hNumerator->Rebin(rebinFactor);
    hDenominator->Rebin(rebinFactor);
    TH1F *hCF = (TH1F*)hNumerator->Clone(name);
    hCF->Divide(hDenominator);
    return hCF;
}


vector<TH1F*> backgroundCorrelation(TFile* filename, TString directoryName, Int_t rebinFactor=1){
    
    TFile* f = filename;
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir = (TDirectory*)f->Get(directoryName);

    TDirectory *sameEvent = (TDirectory*)Dir->Get("SameEvent");
    TDirectory *mixedEvent = (TDirectory*)Dir->Get("MixedEvent");

    TH2F *hAllCentrSameEvent = (TH2F*)sameEvent->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEvent = (TH2F*)mixedEvent->Get("relPairkstarMult");

    TH1D *hSameEventRatio = hAllCentrSameEvent->ProjectionX("hSameEventRatio",2,9);
    TH1D *hMixedEventRatio = hAllCentrMixedEvent->ProjectionX("hMixedEventRatio",2,9);

    TH1D *hkstarSameEvent0_10 = hAllCentrSameEvent->ProjectionX("hkstarSameEvent0_10",2,2);
    TH1D *hkstarMixedEvent0_10 = hAllCentrMixedEvent->ProjectionX("hkstarMixedEvent0_10",2,2);

    TH1D *hkstarSameEvent10_20 = hAllCentrSameEvent->ProjectionX("hkstarSameEvent10_20",3,3);
    TH1D *hkstarMixedEvent10_20 = hAllCentrMixedEvent->ProjectionX("hkstarMixedEvent10_20",3,3);

    TH1D *hkstarSameEvent20_30 = hAllCentrSameEvent->ProjectionX("hkstarSameEvent20_30",4,4);
    TH1D *hkstarMixedEvent20_30 = hAllCentrMixedEvent->ProjectionX("hkstarMixedEvent20_30",4,4);

    TH1D *hkstarSameEvent30_50 = hAllCentrSameEvent->ProjectionX("hkstarSameEvent30_50",5,6);
    TH1D *hkstarMixedEvent30_50 = hAllCentrMixedEvent->ProjectionX("hkstarMixedEvent30_50",5,6);

    TH1D *hkstarSameEvent50_100 = hAllCentrSameEvent->ProjectionX("hkstarSameEvent50_80",7,9);
    TH1D *hkstarMixedEvent50_100 = hAllCentrMixedEvent->ProjectionX("hkstarMixedEvent50_80",7,9);
    

    TH1F* hCFall = createCFHistogram((TH1F*)hSameEventRatio, (TH1F*)hMixedEventRatio, "hCFun", rebinFactor);   //antyproton omega
    TH1F* hCF0_10 = createCFHistogram((TH1F*)hkstarSameEvent0_10, (TH1F*)hkstarMixedEvent0_10, "hCF0_10", rebinFactor);
    TH1F* hCF10_20 = createCFHistogram((TH1F*)hkstarSameEvent10_20, (TH1F*)hkstarMixedEvent10_20, "hCF10_20", rebinFactor);
    TH1F* hCF20_30 = createCFHistogram((TH1F*)hkstarSameEvent20_30, (TH1F*)hkstarMixedEvent20_30, "hCF20_30", rebinFactor);
    TH1F* hCF30_50 = createCFHistogram((TH1F*)hkstarSameEvent30_50, (TH1F*)hkstarMixedEvent30_50, "hCF30_50", rebinFactor); 
    TH1F* hCF50_100 = createCFHistogram((TH1F*)hkstarSameEvent50_100, (TH1F*)hkstarMixedEvent50_100, "hCF50_100", rebinFactor);


    vector<TH1F*> results;
    results.push_back(hCFall);
    results.push_back(hCF0_10);
    results.push_back(hCF10_20);
    results.push_back(hCF20_30);
    results.push_back(hCF30_50);
    results.push_back(hCF50_100);
    results.push_back(hSameEventRatio);
    results.push_back(hMixedEventRatio);
    return results;
}

