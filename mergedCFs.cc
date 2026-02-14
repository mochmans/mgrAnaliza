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
#include "histRead.cc"

void configureHistogram(TH1* histo, Color_t markerColor, Style_t markerStyle, Size_t markerSize){
    histo->SetMarkerColor(markerColor);
    histo->SetMarkerStyle(markerStyle);
    histo->SetMarkerSize(markerSize);
    histo->GetXaxis()->SetRangeUser(0,1.5);
    histo->Scale(1.0);
}

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

HistogramGroups kStarToCFs(HistogramGroups& hists, Int_t rebinFactor=1){
    HistogramGroups cfHists;

    for(size_t i = 0; i < hists.omegaCor.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.omegaCor[i], hists.omegaCor[i+1], Form("cfomega_%zu", i/2), rebinFactor);
        cfHists.omegaCor.push_back(cf);
    }

    for(size_t i = 0; i < hists.antiOmegaCor.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.antiOmegaCor[i], hists.antiOmegaCor[i+1], Form("cf_antyomega_%zu", i/2), rebinFactor);
        cfHists.antiOmegaCor.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundOmegaSum.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundOmegaSum[i], hists.backgroundOmegaSum[i+1], Form("cf_bkg_omega_sum_%zu", i/2), rebinFactor);
        cfHists.backgroundOmegaSum.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundAntiOmegaSum.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundAntiOmegaSum[i], hists.backgroundAntiOmegaSum[i+1], Form("cf_bkg_antyomega_sum_%zu", i/2), rebinFactor);
        cfHists.backgroundAntiOmegaSum.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundHighOmega.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundHighOmega[i], hists.backgroundHighOmega[i+1], Form("cf_bkg_high_omega_%zu", i/2), rebinFactor);
        cfHists.backgroundHighOmega.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundHighAntiOmega.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundHighAntiOmega[i], hists.backgroundHighAntiOmega[i+1], Form("cf_bkg_high_antyomega_%zu", i/2), rebinFactor);
        cfHists.backgroundHighAntiOmega.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundLowOmega.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundLowOmega[i], hists.backgroundLowOmega[i+1], Form("cf_bkg_low_omega_%zu", i/2), rebinFactor);
        cfHists.backgroundLowOmega.push_back(cf);
    }

    for(size_t i = 0; i < hists.backgroundLowAntiOmega.size(); i+=2){
        TH1F* cf = createCFHistogram(hists.backgroundLowAntiOmega[i], hists.backgroundLowAntiOmega[i+1], Form("cf_bkg_low_antyomega_%zu", i/2), rebinFactor);
        cfHists.backgroundLowAntiOmega.push_back(cf);
    }

    return cfHists;
}

void mergedCFs(){


    std::cout<< "Loading histograms from file..." << std::endl;
    HistogramGroups mergedHists = loadkStarHistograms("kStarMergedNozzm.root");

    std::cout<< mergedHists.omegaCor.size() << std::endl; // assuming all have the same number of histograms

    int rebinFactor = 1;

    HistogramGroups mergedCFs = kStarToCFs(mergedHists, rebinFactor);

    std::cout<< "Merged histograms created" << std::endl;

    TCanvas *c1 = new TCanvas("pbarOmega","pbarOmega",800,600);
    c1->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(mergedCFs.omegaCor[i],kBlue,kFullCircle,0.5);
        c1->cd(i+1);
        mergedCFs.omegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFs.omegaCor[i]->Draw();
    }

    TCanvas *c2 = new TCanvas("pOmegabar","pOmegabar",800,600);
    c2->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(mergedCFs.antiOmegaCor[i],kRed,kFullSquare,0.5);
        c2->cd(i+1);
        mergedCFs.antiOmegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFs.antiOmegaCor[i]->Draw();
    }    



    TCanvas *c10 = new TCanvas("Background high","Background high",800,600);
    c10->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(mergedCFs.backgroundHighOmega[i], kBlue, kFullCircle, 0.5);
        configureHistogram(mergedCFs.backgroundHighAntiOmega[i], kRed, kFullSquare, 0.5);
        mergedCFs.backgroundHighOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        c10->cd(i+1);
        mergedCFs.backgroundHighOmega[i]->Draw();
        mergedCFs.backgroundHighAntiOmega[i]->Draw("SAME");
    }
    c10->cd(1);
    TLegend *legendHigh = new TLegend();
    legendHigh->AddEntry(mergedCFs.backgroundHighOmega[0],"Omega background","p");
    legendHigh->AddEntry(mergedCFs.backgroundHighAntiOmega[0],"Antyomega background","p");
    legendHigh->Draw();


    TCanvas *c11 = new TCanvas("Background low","Background low",800,600);
    c11->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(mergedCFs.backgroundLowOmega[i], kBlue, kFullCircle, 0.5);
        configureHistogram(mergedCFs.backgroundLowAntiOmega[i], kRed, kFullSquare, 0.5);
        mergedCFs.backgroundLowOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        c11->cd(i+1);
        mergedCFs.backgroundLowOmega[i]->Draw();
        mergedCFs.backgroundLowAntiOmega[i]->Draw("SAME");
    }
    c11->cd(1);
    TLegend *legendLow = new TLegend();
    legendLow->AddEntry(mergedCFs.backgroundLowOmega[0],"Omega background","p");
    legendLow->AddEntry(mergedCFs.backgroundLowAntiOmega[0],"Antyomega background","p");
    legendLow->Draw();



    TCanvas *c5 = new TCanvas("Cf comparison","Cf comparison",800,600);
    c5->Divide(2,3);
    for(Int_t i = 0; i < 6; i++){
        c5->cd(i+1);
        mergedCFs.omegaCor[i]->Draw();
        mergedCFs.antiOmegaCor[i]->Draw("SAME");
    }

    TLegend *legend3 = new TLegend();
    legend3->AddEntry(mergedCFs.omegaCor[0],"antiproton-omega","p");
    legend3->AddEntry(mergedCFs.antiOmegaCor[0],"proton-antiomega","p");
    c5->cd(1);
    legend3->Draw();



}