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
#include "TDirectory.h"

void configureHistogram(TH1* histo, Color_t markerColor, Style_t markerStyle, Size_t markerSize){
    histo->SetMarkerColor(markerColor);
    histo->SetMarkerStyle(markerStyle);
    histo->SetMarkerSize(markerSize);
    histo->GetXaxis()->SetRangeUser(0,0.4);
    histo->Scale(1.0);
}

TH1F *createCFHistogram(TH1F* hNumerator, TH1F* hDenominator, const char* name){
    TH1F *hCF = (TH1F*)hNumerator->Clone(name);
    hCF->Divide(hDenominator);
    return hCF;
}

void normalizeHistogram(TH1F* histo, Double_t kstarLower, Double_t kstarUpper){
    Int_t binLow = histo->GetXaxis()->FindBin(kstarLower);
    Int_t binUp = histo->GetXaxis()->FindBin(kstarUpper);
    Double_t binAverage = 0;
    for(Int_t i = binLow; i <= binUp; i++){
        binAverage += histo->GetBinContent(i);
    }
    binAverage /= (binUp - binLow + 1);
    histo->Scale(1.0/binAverage);

}

int ratioFull(){
    TH1::SetDefaultSumw2();

    TFile* f = new TFile("mediumDataOmegaAR.root","open");
    TDirectory *pOmegaBar = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_antyomega"); //proton antyomega
    TDirectory *sameEventpAo = (TDirectory*)pOmegaBar->Get("SameEvent");
    TDirectory *mixedEventpAo = (TDirectory*)pOmegaBar->Get("MixedEvent");
    
    

    TH2F *hAllCentrSameEventpAo = (TH2F*)sameEventpAo->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEventpAo = (TH2F*)mixedEventpAo->Get("relPairkstarMult");

    TH1D *hSameEventRatiopAo = hAllCentrSameEventpAo->ProjectionX("hSameEventRatiopAo",2,9);
    TH1D *hMixedEventRatiopAo = hAllCentrMixedEventpAo->ProjectionX("hMixedEventRatiopAo",2,9);

    TH1D *hkstarSameEventpAo0_10 = hAllCentrSameEventpAo->ProjectionX("hkstarSameEventpAo0_10",2,2);
    TH1D *hkstarMixedEventpAo0_10 = hAllCentrMixedEventpAo->ProjectionX("hkstarMixedEventpAo0_10",2,2);

    TH1D *hkstarSameEventpAo10_20 = hAllCentrSameEventpAo->ProjectionX("hkstarSameEventpAo10_20",3,3);
    TH1D *hkstarMixedEventpAo10_20 = hAllCentrMixedEventpAo->ProjectionX("hkstarMixedEventpAo10_20",3,3);   
    TH1D *hkstarSameEventpAo20_30 = hAllCentrSameEventpAo->ProjectionX("hkstarSameEventpAo20_30",4,4);
    TH1D *hkstarMixedEventpAo20_30 = hAllCentrMixedEventpAo->ProjectionX("hkstarMixedEventpAo20_30",4,4);

    TH1D *hkstarSameEventpAo30_50 = hAllCentrSameEventpAo->ProjectionX("hkstarSameEventpAo30_50",5,6);
    TH1D *hkstarMixedEventpAo30_50 = hAllCentrMixedEventpAo->ProjectionX("hkstarMixedEventpAo30_50",5,6);  
    TH1D *hkstarSameEventpAo50_100 = hAllCentrSameEventpAo->ProjectionX("hkstarSameEventpAo50_100",7,9);
    TH1D *hkstarMixedEventpAo50_100 = hAllCentrMixedEventpAo->ProjectionX("hkstarMixedEventpAo50_100",7,9);


    TDirectory *pbarOmega = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended"); //antyproton omega
    TDirectory *sameEvent = (TDirectory*)pbarOmega->Get("SameEvent");
    TDirectory *mixedEvent = (TDirectory*)pbarOmega->Get("MixedEvent");


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



   
    TH1F* hCFall = createCFHistogram((TH1F*)hSameEventRatio, (TH1F*)hMixedEventRatio, "hCFunpbar");   //antyproton omega
    
    TH1F* hCF0_10 = createCFHistogram((TH1F*)hkstarSameEvent0_10, (TH1F*)hkstarMixedEvent0_10, "hCF0_10pbar");
    
    TH1F* hCF10_20 = createCFHistogram((TH1F*)hkstarSameEvent10_20, (TH1F*)hkstarMixedEvent10_20, "hCF10_20pbar");
    
    TH1F* hCF20_30 = createCFHistogram((TH1F*)hkstarSameEvent20_30, (TH1F*)hkstarMixedEvent20_30, "hCF20_30pbar");
    
    TH1F* hCF30_50 = createCFHistogram((TH1F*)hkstarSameEvent30_50, (TH1F*)hkstarMixedEvent30_50, "hCF30_50pbar");
    
    TH1F* hCF50_100 = createCFHistogram((TH1F*)hkstarSameEvent50_100, (TH1F*)hkstarMixedEvent50_100, "hCF50_100pbar");
    


    TH1F* hCFallpAo = createCFHistogram((TH1F*)hSameEventRatiopAo, (TH1F*)hMixedEventRatiopAo, "hCFunp");  //proton Antyomega
    
    TH1F* hCF0_10pAo = createCFHistogram((TH1F*)hkstarSameEventpAo0_10, (TH1F*)hkstarMixedEventpAo0_10, "hCF0_10p");
    
    TH1F* hCF10_20pAo = createCFHistogram((TH1F*)hkstarSameEventpAo10_20, (TH1F*)hkstarMixedEventpAo10_20, "hCF10_20p");
    
    TH1F* hCF20_30pAo = createCFHistogram((TH1F*)hkstarSameEventpAo20_30, (TH1F*)hkstarMixedEventpAo20_30, "hCF20_30p");
    
    TH1F* hCF30_50pAo = createCFHistogram((TH1F*)hkstarSameEventpAo30_50, (TH1F*)hkstarMixedEventpAo30_50, "hCF30_50p");
    
    TH1F* hCF50_100pAo = createCFHistogram((TH1F*)hkstarSameEventpAo50_100, (TH1F*)hkstarMixedEventpAo50_100, "hCF50_100p");
     

    Int_t rebinF = 4;
    hCFall->Rebin(rebinF);
    hCF0_10->Rebin(rebinF);
    hCF10_20->Rebin(rebinF);
    hCF20_30->Rebin(rebinF);
    hCF30_50->Rebin(rebinF);
    hCF50_100->Rebin(rebinF);
    hCFallpAo->Rebin(rebinF);
    hCF0_10pAo->Rebin(rebinF);
    hCF10_20pAo->Rebin(rebinF);
    hCF20_30pAo->Rebin(rebinF);
    hCF30_50pAo->Rebin(rebinF);
    hCF50_100pAo->Rebin(rebinF);

    configureHistogram(hCFall, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF0_10, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF10_20, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF20_30, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF30_50, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF50_100, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCFallpAo, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF0_10pAo, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF10_20pAo, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF20_30pAo, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF30_50pAo, kBlue, kOpenCircle, 0.5);
    configureHistogram(hCF50_100pAo, kBlue, kOpenCircle, 0.5);

    normalizeHistogram(hCFall, 0.5, 1.5);
    normalizeHistogram(hCF0_10, 0.5, 1.5);
    normalizeHistogram(hCF10_20, 0.5, 1.5);
    normalizeHistogram(hCF20_30, 0.5, 1.5); 
    normalizeHistogram(hCF30_50, 0.5, 1.5);
    normalizeHistogram(hCF50_100, 0.5, 1.5);
    normalizeHistogram(hCFallpAo, 0.5, 1.5);
    normalizeHistogram(hCF0_10pAo, 0.5, 1.5);
    normalizeHistogram(hCF10_20pAo, 0.5, 1.5);
    normalizeHistogram(hCF20_30pAo, 0.5, 1.5); 
    normalizeHistogram(hCF30_50pAo, 0.5, 1.5);
    normalizeHistogram(hCF50_100pAo, 0.5, 1.5);

    TCanvas *c1 = new TCanvas("pbarOmega","pbarOmega",800,600);
    c1->Divide(2,3);
    c1->cd(1);
    hCFall->Draw();
    c1->cd(2);
    hCF0_10->Draw();
    c1->cd(3);
    hCF10_20->Draw();
    c1->cd(4);
    hCF20_30->Draw();
    c1->cd(5);
    hCF30_50->Draw();
    c1->cd(6);
    hCF50_100->Draw();

    TCanvas *c2 = new TCanvas("pOmegabar","pOmegabar",800,600);
    c2->Divide(2,3);
    c2->cd(1);
    hCFallpAo->Draw();
    c2->cd(2);
    hCF0_10pAo->Draw();
    c2->cd(3);
    hCF10_20pAo->Draw();
    c2->cd(4);
    hCF20_30pAo->Draw();
    c2->cd(5);
    hCF30_50pAo->Draw();
    c2->cd(6);
    hCF50_100pAo->Draw(); 

    
    return 0; 
}