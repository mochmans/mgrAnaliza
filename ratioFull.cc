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

void configureHistogram(TH1* histo, Color_t markerColor, Style_t markerStyle, Size_t markerSize){
    histo->SetMarkerColor(markerColor);
    histo->SetMarkerStyle(markerStyle);
    histo->SetMarkerSize(markerSize);
    histo->GetXaxis()->SetRangeUser(0,1.5);
    histo->Scale(1.0);
}

TH1F *createCFHistogram(TH1F* hNumerator, TH1F* hDenominator, const char* name, Int_t rebinFactor=1){
    // hNumerator->Rebin(rebinFactor);
    // hDenominator->Rebin(rebinFactor);
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

    TFile* f = new TFile("mediumDataOmegaARfull.root","open");
    


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


    Int_t rebinFactor = 3;

    hkstarMixedEvent0_10->Rebin(rebinFactor);
    hkstarMixedEvent10_20->Rebin(rebinFactor);
    hkstarMixedEvent20_30->Rebin(rebinFactor);
    hkstarMixedEvent30_50->Rebin(rebinFactor);
    hkstarMixedEvent50_100->Rebin(rebinFactor);
    hMixedEventRatio->Rebin(rebinFactor);
    hkstarMixedEventpAo0_10->Rebin(rebinFactor);
    hkstarMixedEventpAo10_20->Rebin(rebinFactor);
    hkstarMixedEventpAo20_30->Rebin(rebinFactor);
    hkstarMixedEventpAo30_50->Rebin(rebinFactor);
    hkstarMixedEventpAo50_100->Rebin(rebinFactor);
    hMixedEventRatiopAo->Rebin(rebinFactor);
    hSameEventRatio->Rebin(rebinFactor);
    hkstarSameEvent0_10->Rebin(rebinFactor);
    hkstarSameEvent10_20->Rebin(rebinFactor);
    hkstarSameEvent20_30->Rebin(rebinFactor);
    hkstarSameEvent30_50->Rebin(rebinFactor);
    hkstarSameEvent50_100->Rebin(rebinFactor);
    hkstarSameEventpAo0_10->Rebin(rebinFactor);
    hkstarSameEventpAo10_20->Rebin(rebinFactor);
    hkstarSameEventpAo20_30->Rebin(rebinFactor);
    hkstarSameEventpAo30_50->Rebin(rebinFactor);
    hkstarSameEventpAo50_100->Rebin(rebinFactor);
    hSameEventRatiopAo->Rebin(rebinFactor);
   
    TH1F* hCFall = createCFHistogram((TH1F*)hSameEventRatio, (TH1F*)hMixedEventRatio, "hCFunpbar", rebinFactor);   //antyproton omega
    
    TH1F* hCF0_10 = createCFHistogram((TH1F*)hkstarSameEvent0_10, (TH1F*)hkstarMixedEvent0_10, "hCF0_10pbar", rebinFactor);
    
    TH1F* hCF10_20 = createCFHistogram((TH1F*)hkstarSameEvent10_20, (TH1F*)hkstarMixedEvent10_20, "hCF10_20pbar", rebinFactor);
    
    TH1F* hCF20_30 = createCFHistogram((TH1F*)hkstarSameEvent20_30, (TH1F*)hkstarMixedEvent20_30, "hCF20_30pbar", rebinFactor);
    
    TH1F* hCF30_50 = createCFHistogram((TH1F*)hkstarSameEvent30_50, (TH1F*)hkstarMixedEvent30_50, "hCF30_50pbar", rebinFactor);
    
    TH1F* hCF50_100 = createCFHistogram((TH1F*)hkstarSameEvent50_100, (TH1F*)hkstarMixedEvent50_100, "hCF50_100pbar", rebinFactor);
    


    TH1F* hCFallpAo = createCFHistogram((TH1F*)hSameEventRatiopAo, (TH1F*)hMixedEventRatiopAo, "hCFunp",rebinFactor);  //proton Antyomega
    
    TH1F* hCF0_10pAo = createCFHistogram((TH1F*)hkstarSameEventpAo0_10, (TH1F*)hkstarMixedEventpAo0_10, "hCF0_10p", rebinFactor);
    
    TH1F* hCF10_20pAo = createCFHistogram((TH1F*)hkstarSameEventpAo10_20, (TH1F*)hkstarMixedEventpAo10_20, "hCF10_20p", rebinFactor);
    
    TH1F* hCF20_30pAo = createCFHistogram((TH1F*)hkstarSameEventpAo20_30, (TH1F*)hkstarMixedEventpAo20_30, "hCF20_30p", rebinFactor);
    
    TH1F* hCF30_50pAo = createCFHistogram((TH1F*)hkstarSameEventpAo30_50, (TH1F*)hkstarMixedEventpAo30_50, "hCF30_50p", rebinFactor);
    
    TH1F* hCF50_100pAo = createCFHistogram((TH1F*)hkstarSameEventpAo50_100, (TH1F*)hkstarMixedEventpAo50_100, "hCF50_100p", rebinFactor);
     


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
    hCFall->GetYaxis()->SetRangeUser(0.4,1.6);
    hCFall->Draw();
    c1->cd(2);
    hCF0_10->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF0_10->Draw();
    c1->cd(3);
    hCF10_20->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF10_20->Draw();
    c1->cd(4);
    hCF20_30->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF20_30->Draw();
    c1->cd(5);
    hCF30_50->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF30_50->Draw();
    c1->cd(6);
    hCF50_100->GetYaxis()->SetRangeUser(0,3);
    hCF50_100->Draw();

    TCanvas *c2 = new TCanvas("pOmegabar","pOmegabar",800,600);
    c2->Divide(2,3);
    c2->cd(1);
    hCFallpAo->GetYaxis()->SetRangeUser(0.4,1.6);
    hCFallpAo->Draw();
    c2->cd(2);
    hCF0_10pAo->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF0_10pAo->Draw();
    c2->cd(3);
    hCF10_20pAo->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF10_20pAo->Draw();
    c2->cd(4);
    hCF20_30pAo->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF20_30pAo->Draw();
    c2->cd(5);
    hCF30_50pAo->GetYaxis()->SetRangeUser(0.4,1.6);
    hCF30_50pAo->Draw();
    c2->cd(6);
    hCF50_100pAo->GetYaxis()->SetRangeUser(0,3);
    hCF50_100pAo->Draw(); 

    double_t msize = 0.25;

    TCanvas *c3 = new TCanvas("events proton antiomega","Comparison",800,600);
    c3->Divide(2,3);
    c3->cd(1);
    hSameEventRatiopAo->SetMarkerColor(kBlue);
    hSameEventRatiopAo->GetXaxis()->SetRangeUser(0,4);
    hSameEventRatiopAo->Draw();
    hMixedEventRatiopAo->SetMarkerStyle(kOpenSquare);
    hMixedEventRatiopAo->SetMarkerColor(kRed);
    hMixedEventRatiopAo->SetMarkerSize(msize);
    hMixedEventRatiopAo->Draw("SAME");

    TLegend *legend1 = new TLegend();
    legend1->AddEntry(hkstarSameEventpAo0_10,"Same Event","p");
    legend1->AddEntry(hkstarMixedEventpAo0_10,"Mixed Event","p");
    legend1->Draw();

    

    c3->cd(2);  
    hkstarSameEventpAo0_10->SetMarkerColor(kBlue);
    hkstarSameEventpAo0_10->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEventpAo0_10->Draw();
    hkstarMixedEventpAo0_10->SetMarkerStyle(kOpenSquare);
    hkstarMixedEventpAo0_10->SetMarkerColor(kRed);
    hkstarMixedEventpAo0_10->SetMarkerSize(msize);
    hkstarMixedEventpAo0_10->Draw("SAME");
    c3->cd(3);
    hkstarSameEventpAo10_20->SetMarkerColor(kBlue);
    hkstarSameEventpAo10_20->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEventpAo10_20->Draw();
    hkstarMixedEventpAo10_20->SetMarkerStyle(kOpenSquare);
    hkstarMixedEventpAo10_20->SetMarkerColor(kRed);
    hkstarMixedEventpAo10_20->SetMarkerSize(msize);
    hkstarMixedEventpAo10_20->Draw("SAME");
    c3->cd(4);
    hkstarSameEventpAo20_30->SetMarkerColor(kBlue);
    hkstarSameEventpAo20_30->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEventpAo20_30->Draw();
    hkstarMixedEventpAo20_30->SetMarkerStyle(kOpenSquare);
    hkstarMixedEventpAo20_30->SetMarkerColor(kRed);
    hkstarMixedEventpAo20_30->SetMarkerSize(msize);
    hkstarMixedEventpAo20_30->Draw("SAME");
    c3->cd(5);
    hkstarSameEventpAo30_50->SetMarkerColor(kBlue);
    hkstarSameEventpAo30_50->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEventpAo30_50->Draw();
    hkstarMixedEventpAo30_50->SetMarkerStyle(kOpenSquare);
    hkstarMixedEventpAo30_50->SetMarkerColor(kRed);
    hkstarMixedEventpAo30_50->SetMarkerSize(msize);
    hkstarMixedEventpAo30_50->Draw("SAME");
    c3->cd(6);
    hkstarSameEventpAo50_100->SetMarkerColor(kBlue);
    hkstarSameEventpAo50_100->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEventpAo50_100->Draw();
    hkstarMixedEventpAo50_100->SetMarkerStyle(kOpenSquare);
    hkstarMixedEventpAo50_100->SetMarkerColor(kRed);
    hkstarMixedEventpAo50_100->SetMarkerSize(msize);
    hkstarMixedEventpAo50_100->Draw("SAME");

    

    TCanvas *c4 = new TCanvas("events antiproton omega","k* Tcomparison",800,600);
    c4->Divide(2,3);
    c4->cd(1);
    hSameEventRatio->SetMarkerColor(kBlue);
    hSameEventRatio->GetXaxis()->SetRangeUser(0,4);
    hSameEventRatio->Draw();
    hMixedEventRatio->SetMarkerStyle(kOpenSquare);
    hMixedEventRatio->SetMarkerColor(kRed);
    hMixedEventRatio->SetMarkerSize(msize);
    hMixedEventRatio->Draw("SAME");
    c4->cd(2);  
    hkstarSameEvent0_10->SetMarkerColor(kBlue);   
    hkstarSameEvent0_10->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEvent0_10->Draw();
    hkstarMixedEvent0_10->SetMarkerStyle(kOpenSquare);
    hkstarMixedEvent0_10->SetMarkerColor(kRed);
    hkstarMixedEvent0_10->SetMarkerSize(msize);
    hkstarMixedEvent0_10->Draw("SAME"); 
    c4->cd(3);
    hkstarSameEvent10_20->SetMarkerColor(kBlue);
    hkstarSameEvent10_20->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEvent10_20->Draw();
    hkstarMixedEvent10_20->SetMarkerStyle(kOpenSquare);
    hkstarMixedEvent10_20->SetMarkerColor(kRed);
    hkstarMixedEvent10_20->SetMarkerSize(msize);
    hkstarMixedEvent10_20->Draw("SAME");
    c4->cd(4);
    hkstarSameEvent20_30->SetMarkerColor(kBlue);
    hkstarSameEvent20_30->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEvent20_30->Draw();
    hkstarMixedEvent20_30->SetMarkerStyle(kOpenSquare);
    hkstarMixedEvent20_30->SetMarkerColor(kRed);
    hkstarMixedEvent20_30->SetMarkerSize(msize);
    hkstarMixedEvent20_30->Draw("SAME");    
    c4->cd(5);
    hkstarSameEvent30_50->SetMarkerColor(kBlue);
    hkstarSameEvent30_50->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEvent30_50->Draw();
    hkstarMixedEvent30_50->SetMarkerStyle(kOpenSquare);
    hkstarMixedEvent30_50->SetMarkerColor(kRed);
    hkstarMixedEvent30_50->SetMarkerSize(msize);
    hkstarMixedEvent30_50->Draw("SAME");
    c4->cd(6);
    hkstarSameEvent50_100->SetMarkerColor(kBlue);
    hkstarSameEvent50_100->GetXaxis()->SetRangeUser(0,4);
    hkstarSameEvent50_100->Draw();
    hkstarMixedEvent50_100->SetMarkerStyle(kOpenSquare);
    hkstarMixedEvent50_100->SetMarkerColor(kRed);
    hkstarMixedEvent50_100->SetMarkerSize(msize);
    hkstarMixedEvent50_100->Draw("SAME");

    TLegend *legend2 = new TLegend();
    legend2->AddEntry(hkstarSameEvent0_10,"Same Event","p");
    legend2->AddEntry(hkstarMixedEvent0_10,"Mixed Event","p");
    c4->cd(1);
    legend2->Draw();

    TCanvas *c5 = new TCanvas("Cf comparison","Cf comparison",800,600);
    c5->Divide(2,3);
    c5->cd(1);
    hCFall->SetMarkerColor(kBlue);
    hCFall->Draw();
    hCFallpAo->SetMarkerStyle(kOpenSquare);
    hCFallpAo->SetMarkerColor(kRed);
    hCFallpAo->Draw("SAME");
    c5->cd(2);
    hCF0_10->SetMarkerColor(kBlue);
    hCF0_10->Draw();
    hCF0_10pAo->SetMarkerStyle(kOpenSquare);
    hCF0_10pAo->SetMarkerColor(kRed);
    hCF0_10pAo->Draw("SAME");
    c5->cd(3);
    hCF10_20->SetMarkerColor(kBlue);
    hCF10_20->Draw();
    hCF10_20pAo->SetMarkerStyle(kOpenSquare);
    hCF10_20pAo->SetMarkerColor(kRed);
    hCF10_20pAo->Draw("SAME");
    c5->cd(4);
    hCF20_30->SetMarkerColor(kBlue);
    hCF20_30->Draw();
    hCF20_30pAo->SetMarkerStyle(kOpenSquare);
    hCF20_30pAo->SetMarkerColor(kRed);
    hCF20_30pAo->Draw("SAME");
    c5->cd(5);
    hCF30_50->SetMarkerColor(kBlue);
    hCF30_50->Draw();
    hCF30_50pAo->SetMarkerStyle(kOpenSquare);
    hCF30_50pAo->SetMarkerColor(kRed);
    hCF30_50pAo->Draw("SAME");
    c5->cd(6);
    hCF50_100->SetMarkerColor(kBlue);
    hCF50_100->Draw();
    hCF50_100pAo->SetMarkerStyle(kOpenSquare);
    hCF50_100pAo->SetMarkerColor(kRed); 
    hCF50_100pAo->Draw("SAME");

    TLegend *legend3 = new TLegend();
    legend3->AddEntry(hCFall,"antiproton-omega","p");
    legend3->AddEntry(hCFallpAo,"proton-antiomega","p");
    c5->cd(1);
    legend3->Draw();


    
    return 0; 
}