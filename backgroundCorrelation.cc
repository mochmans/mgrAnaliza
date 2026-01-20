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

void normalizeHistogram(TH1F* histo, Double_t kstarLower=0.5, Double_t kstarUpper=1.5){
    Int_t binLow = histo->GetXaxis()->FindBin(kstarLower);
    Int_t binUp = histo->GetXaxis()->FindBin(kstarUpper);
    Double_t binAverage = 0;
    for(Int_t i = binLow; i <= binUp; i++){
        binAverage += histo->GetBinContent(i);
    }
    binAverage /= (binUp - binLow + 1);
    histo->Scale(1.0/binAverage);

}

TH1F *createCFHistogram(TH1F* hNumerator, TH1F* hDenominator, const char* name, Int_t rebinFactor=1){
    hNumerator->Rebin(rebinFactor);
    hDenominator->Rebin(rebinFactor);
    TH1F *hCF = (TH1F*)hNumerator->Clone(name);
    hCF->Divide(hDenominator);
    normalizeHistogram(hCF,0.5,1.5);
    return hCF;
}


std::vector<TH1F*> correlationFunctions(TFile* filename, TString directoryName, Int_t rebinFactor=1){
    
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


    std::vector<TH1F*> results;
    results.push_back(hCFall);
    results.push_back(hCF0_10);
    results.push_back(hCF10_20);
    results.push_back(hCF20_30);
    results.push_back(hCF30_50);
    results.push_back(hCF50_100);
    results.push_back((TH1F*)hSameEventRatio);
    results.push_back((TH1F*)hMixedEventRatio);
    return results;
}

std::vector<TH1F*> backgroundCorrelationSum(TFile* filename, TString directoryNameLow ,TString directoryNameHigh, Int_t rebinFactor=1){
    
    TFile* f = filename;
    
    TH1::SetDefaultSumw2();

    TDirectory *DirLow = (TDirectory*)f->Get(directoryNameLow);
    TDirectory *DirHigh = (TDirectory*)f->Get(directoryNameHigh);

    TDirectory *sameEventLow = (TDirectory*)DirLow->Get("SameEvent");
    TDirectory *mixedEventLow = (TDirectory*)DirLow->Get("MixedEvent");
    TDirectory *sameEventHigh = (TDirectory*)DirHigh->Get("SameEvent");
    TDirectory *mixedEventHigh = (TDirectory*)DirHigh->Get("MixedEvent");

    TH2F *hAllCentrSameEventLow = (TH2F*)sameEventLow->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEventLow = (TH2F*)mixedEventLow->Get("relPairkstarMult");
    TH2F *hAllCentrSameEventHigh = (TH2F*)sameEventHigh->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEventHigh = (TH2F*)mixedEventHigh->Get("relPairkstarMult");


    TH1D *hSameEventRatioLow = hAllCentrSameEventLow->ProjectionX("hSameEventRatioLow",2,9);
    TH1D *hMixedEventRatioLow = hAllCentrMixedEventLow->ProjectionX("hMixedEventRatioLow",2,9);

    TH1D *hkstarSameEventLow0_10 = hAllCentrSameEventLow->ProjectionX("hkstarSameEventLow0_10",2,2);
    TH1D *hkstarMixedEventLow0_10 = hAllCentrMixedEventLow->ProjectionX("hkstarMixedEventLow0_10",2,2);   
    TH1D *hkstarSameEventLow10_20 = hAllCentrSameEventLow->ProjectionX("hkstarSameEventLow10_20",3,3);
    TH1D *hkstarMixedEventLow10_20 = hAllCentrMixedEventLow->ProjectionX("hkstarMixedEventLow10_20",3,3);

    TH1D *hkstarSameEventLow20_30 = hAllCentrSameEventLow->ProjectionX("hkstarSameEventLow20_30",4,4);
    TH1D *hkstarMixedEventLow20_30 = hAllCentrMixedEventLow->ProjectionX("hkstarMixedEventLow20_30",4,4);

    TH1D *hkstarSameEventLow30_50 = hAllCentrSameEventLow->ProjectionX("hkstarSameEventLow30_50",5,6);
    TH1D *hkstarMixedEventLow30_50 = hAllCentrMixedEventLow->ProjectionX("hkstarMixedEventLow30_50",5,6);

    TH1D *hkstarSameEventLow50_100 = hAllCentrSameEventLow->ProjectionX("hkstarSameEventLow50_100",7,9);
    TH1D *hkstarMixedEventLow50_100 = hAllCentrMixedEventLow->ProjectionX("hkstarMixedEventLow50_100",7,9);

    TH1D *hSameEventRatioHigh = hAllCentrSameEventHigh->ProjectionX("hSameEventRatioHigh",2,9);
    TH1D *hMixedEventRatioHigh = hAllCentrMixedEventHigh->ProjectionX("hMixedEventRatioHigh",2,9);  
    TH1D *hkstarSameEventHigh0_10 = hAllCentrSameEventHigh->ProjectionX("hkstarSameEventHigh0_10",2,2);
    TH1D *hkstarMixedEventHigh0_10 = hAllCentrMixedEventHigh->ProjectionX("hkstarMixedEventHigh0_10",2,2);
    TH1D *hkstarSameEventHigh10_20 = hAllCentrSameEventHigh->ProjectionX("hkstarSameEventHigh10_20",3,3);
    TH1D *hkstarMixedEventHigh10_20 = hAllCentrMixedEventHigh->ProjectionX("hkstarMixedEventHigh10_20",3,3);
    TH1D *hkstarSameEventHigh20_30 = hAllCentrSameEventHigh->ProjectionX("hkstarSameEventHigh20_30",4,4);
    TH1D *hkstarMixedEventHigh20_30 = hAllCentrMixedEventHigh->ProjectionX("hkstarMixedEventHigh20_30",4,4);
    TH1D *hkstarSameEventHigh30_50 = hAllCentrSameEventHigh->ProjectionX("hkstarSameEventHigh30_50",5,6);
    TH1D *hkstarMixedEventHigh30_50 = hAllCentrMixedEventHigh->ProjectionX("hkstarMixedEventHigh30_50",5,6);
    TH1D *hkstarSameEventHigh50_100 = hAllCentrSameEventHigh->ProjectionX("hkstarSameEventHigh50_100",7,9);
    TH1D *hkstarMixedEventHigh50_100 = hAllCentrMixedEventHigh->ProjectionX("hkstarMixedEventHigh50_100",7,9);

    hSameEventRatioLow->Add(hSameEventRatioHigh);
    hMixedEventRatioLow->Add(hMixedEventRatioHigh);
    hkstarSameEventLow0_10->Add(hkstarSameEventHigh0_10);
    hkstarMixedEventLow0_10->Add(hkstarMixedEventHigh0_10);
    hkstarSameEventLow10_20->Add(hkstarSameEventHigh10_20);
    hkstarMixedEventLow10_20->Add(hkstarMixedEventHigh10_20);
    hkstarSameEventLow20_30->Add(hkstarSameEventHigh20_30);
    hkstarMixedEventLow20_30->Add(hkstarMixedEventHigh20_30);
    hkstarSameEventLow30_50->Add(hkstarSameEventHigh30_50);
    hkstarMixedEventLow30_50->Add(hkstarMixedEventHigh30_50);
    hkstarSameEventLow50_100->Add(hkstarSameEventHigh50_100);
    hkstarMixedEventLow50_100->Add(hkstarMixedEventHigh50_100);
    

    TH1F* hCFall = createCFHistogram((TH1F*)hSameEventRatioLow, (TH1F*)hMixedEventRatioLow, "hCFall", rebinFactor);   //antyproton omega
    TH1F* hCF0_10 = createCFHistogram((TH1F*)hkstarSameEventLow0_10, (TH1F*)hkstarMixedEventLow0_10, "hCF0_10", rebinFactor);
    TH1F* hCF10_20 = createCFHistogram((TH1F*)hkstarSameEventLow10_20, (TH1F*)hkstarMixedEventLow10_20, "hCF10_20", rebinFactor);
    TH1F* hCF20_30 = createCFHistogram((TH1F*)hkstarSameEventLow20_30, (TH1F*)hkstarMixedEventLow20_30, "hCF20_30", rebinFactor);
    TH1F* hCF30_50 = createCFHistogram((TH1F*)hkstarSameEventLow30_50, (TH1F*)hkstarMixedEventLow30_50, "hCF30_50", rebinFactor); 
    TH1F* hCF50_100 = createCFHistogram((TH1F*)hkstarSameEventLow50_100, (TH1F*)hkstarMixedEventLow50_100, "hCF50_100", rebinFactor);


    std::vector<TH1F*> results;
    results.push_back(hCFall);
    results.push_back(hCF0_10);
    results.push_back(hCF10_20);
    results.push_back(hCF20_30);
    results.push_back(hCF30_50);
    results.push_back(hCF50_100);
    results.push_back((TH1F*)hSameEventRatioLow);
    results.push_back((TH1F*)hMixedEventRatioLow);
    return results;
}

std::vector<TH1F*> kstarDistributions(TFile* filename, TString directoryName, Int_t rebinFactor=1){
    
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

    std::vector<TH1F*> results;
    results.push_back((TH1F*)hSameEventRatio);
    results.push_back((TH1F*)hMixedEventRatio);
    results.push_back((TH1F*)hkstarSameEvent0_10);
    results.push_back((TH1F*)hkstarMixedEvent0_10);
    results.push_back((TH1F*)hkstarSameEvent10_20);
    results.push_back((TH1F*)hkstarMixedEvent10_20);
    results.push_back((TH1F*)hkstarSameEvent20_30);
    results.push_back((TH1F*)hkstarMixedEvent20_30);
    results.push_back((TH1F*)hkstarSameEvent30_50);
    results.push_back((TH1F*)hkstarMixedEvent30_50);
    results.push_back((TH1F*)hkstarSameEvent50_100);
    results.push_back((TH1F*)hkstarMixedEvent50_100);
    return results;
}