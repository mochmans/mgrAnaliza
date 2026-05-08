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
    histo->SetLineColor(markerColor);
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

HistogramGroups kStarToCFs(HistogramGroups hists, Int_t rebinFactor=1){
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

HistogramGroups kStarMergedPairAntipair(HistogramGroups hists){
 HistogramGroups mergedHists;

    
    mergedHists.omegaCor.push_back((TH1F*)hists.omegaCor[0]->Clone("mergedOmegakStar"));
    mergedHists.omegaCor[0]->Add(hists.antiOmegaCor[0]);
        for(int j = 1; j < hists.omegaCor.size(); j++) {
            TH1F* temp = (TH1F*)hists.omegaCor[j]->Clone(Form("mergedOmegakStar_%d", j));
            temp->Add(hists.antiOmegaCor[j]);
            mergedHists.omegaCor.push_back(temp);
        }
    

    mergedHists.backgroundLowOmega.push_back((TH1F*)hists.backgroundLowOmega[0]->Clone("mergedLowBgOmegakStar"));
    mergedHists.backgroundLowOmega[0]->Add(hists.backgroundLowAntiOmega[0]);
        for(int j = 1; j < hists.backgroundLowOmega.size(); j++) {
            TH1F* temp = (TH1F*)hists.backgroundLowOmega[j]->Clone(Form("mergedLowBgOmegakStar_%d", j));
            temp->Add(hists.backgroundLowAntiOmega[j]);
            mergedHists.backgroundLowOmega.push_back(temp);
        }
    mergedHists.backgroundHighOmega.push_back((TH1F*)hists.backgroundHighOmega[0]->Clone("mergedHighBgOmegakStar"));
    mergedHists.backgroundHighOmega[0]->Add(hists.backgroundHighAntiOmega[0]);
        for(int j = 1; j < hists.backgroundHighOmega.size(); j++) {
            TH1F* temp = (TH1F*)hists.backgroundHighOmega[j]->Clone(Form("mergedHighBgOmegakStar_%d", j));
            temp->Add(hists.backgroundHighAntiOmega[j]);
            mergedHists.backgroundHighOmega.push_back(temp);
        }
    mergedHists.backgroundOmegaSum.push_back((TH1F*)hists.backgroundOmegaSum[0]->Clone("mergedSumBgOmegakStar"));
    mergedHists.backgroundOmegaSum[0]->Add(hists.backgroundAntiOmegaSum[0]);
        for(int j = 1; j < hists.backgroundOmegaSum.size(); j++) {
            TH1F* temp = (TH1F*)hists.backgroundOmegaSum[j]->Clone(Form("mergedSumBgOmegakStar_%d", j));
            temp->Add(hists.backgroundAntiOmegaSum[j]);
            mergedHists.backgroundOmegaSum.push_back(temp);
        }

    return mergedHists;
}

HistogramGroups kStarToCFsMergedpairAntipair(HistogramGroups hists, Int_t rebinFactor=1){
    HistogramGroups cfHists;

    std::vector<TH1F*> mergedSameEventVec;
    mergedSameEventVec.push_back((TH1F*)hists.omegaCor[0]->Clone("mergedOmegakStar"));
    mergedSameEventVec[0]->Add(hists.antiOmegaCor[0]);
        for(int j = 2; j < hists.omegaCor.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.omegaCor[j]->Clone(Form("mergedOmegakStar_%d", j));
            temp->Add(hists.antiOmegaCor[j]);
            mergedSameEventVec.push_back(temp);
        }
    
    std::vector<TH1F*> mergedMixedEventVec;
    mergedMixedEventVec.push_back((TH1F*)hists.omegaCor[1]->Clone("mergedOmegakStarMixed"));
    mergedMixedEventVec[0]->Add(hists.antiOmegaCor[1]);
        for(int j = 3; j < hists.omegaCor.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.omegaCor[j]->Clone(Form("mergedOmegakStarMixed_%d", j));
            temp->Add(hists.antiOmegaCor[j]);
            mergedMixedEventVec.push_back(temp);
        }
    
    for(size_t i = 0; i < hists.omegaCor.size()/2; i++){
        TH1F* cf = createCFHistogram(mergedSameEventVec[i], mergedMixedEventVec[i], Form("cf_merged_pair_antipair_%zu", i), rebinFactor);
        cfHists.omegaCor.push_back(cf);
    }

    std::vector<TH1F*> mergedSameEventBkgHighVec;
    mergedSameEventBkgHighVec.push_back((TH1F*)hists.backgroundHighOmega[0]->Clone("mergedHighBkgOmega"));
    mergedSameEventBkgHighVec[0]->Add(hists.backgroundHighAntiOmega[0]);
        for(int j = 2; j < hists.backgroundHighOmega.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundHighOmega[j]->Clone(Form("mergedHighBkgOmega_%d", j));
            temp->Add(hists.backgroundHighAntiOmega[j]);
            mergedSameEventBkgHighVec.push_back(temp);
        }
    std::vector<TH1F*> mergedMixedEventBkgHighVec;
    mergedMixedEventBkgHighVec.push_back((TH1F*)hists.backgroundHighOmega[1]->Clone("mergedHighBkgOmegaMixed"));
    mergedMixedEventBkgHighVec[0]->Add(hists.backgroundHighAntiOmega[1]);
        for(int j = 3; j < hists.backgroundHighOmega.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundHighOmega[j]->Clone(Form("mergedHighBkgOmegaMixed_%d", j));
            temp->Add(hists.backgroundHighAntiOmega[j]);
            mergedMixedEventBkgHighVec.push_back(temp);
        }
    for(size_t i = 0; i < hists.backgroundHighOmega.size()/2; i++){
        TH1F* cf = createCFHistogram(mergedSameEventBkgHighVec[i], mergedMixedEventBkgHighVec[i], Form("cf_merged_high_bkg_pair_antipair_%zu", i), rebinFactor);
        cfHists.backgroundHighOmega.push_back(cf);
    }


    std::vector<TH1F*> mergedSameEventBkgLowVec;
    mergedSameEventBkgLowVec.push_back((TH1F*)hists.backgroundLowOmega[0]->Clone("mergedLowBkgOmega"));
    mergedSameEventBkgLowVec[0]->Add(hists.backgroundLowAntiOmega[0]);
        for(int j = 2; j < hists.backgroundLowOmega.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundLowOmega[j]->Clone(Form("mergedLowBkgOmega_%d", j));
            temp->Add(hists.backgroundLowAntiOmega[j]);
            mergedSameEventBkgLowVec.push_back(temp);
        }
    std::vector<TH1F*> mergedMixedEventBkgLowVec;
    mergedMixedEventBkgLowVec.push_back((TH1F*)hists.backgroundLowOmega[1]->Clone("mergedLowBkgOmegaMixed"));
    mergedMixedEventBkgLowVec[0]->Add(hists.backgroundLowAntiOmega[1]);
        for(int j = 3; j < hists.backgroundLowOmega.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundLowOmega[j]->Clone(Form("mergedLowBkgOmegaMixed_%d", j));
            temp->Add(hists.backgroundLowAntiOmega[j]);
            mergedMixedEventBkgLowVec.push_back(temp);
        }
    for(size_t i = 0; i < hists.backgroundLowOmega.size()/2; i++){
        TH1F* cf = createCFHistogram(mergedSameEventBkgLowVec[i], mergedMixedEventBkgLowVec[i], Form("cf_merged_low_bkg_pair_antipair_%zu", i), rebinFactor);
        cfHists.backgroundLowOmega.push_back(cf);
    }

    std::vector<TH1F*> mergedSameEventBkgSumVec;
    mergedSameEventBkgSumVec.push_back((TH1F*)hists.backgroundOmegaSum[0]->Clone("mergedBkgSumOmega"));
    mergedSameEventBkgSumVec[0]->Add(hists.backgroundAntiOmegaSum[0]);
        for(int j = 2; j < hists.backgroundOmegaSum.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundOmegaSum[j]->Clone(Form("mergedBkgSumOmega_%d", j));
            temp->Add(hists.backgroundAntiOmegaSum[j]);
            mergedSameEventBkgSumVec.push_back(temp);
        }
    std::vector<TH1F*> mergedMixedEventBkgSumVec;
    mergedMixedEventBkgSumVec.push_back((TH1F*)hists.backgroundOmegaSum[1]->Clone("mergedBkgSumOmegaMixed"));
    mergedMixedEventBkgSumVec[0]->Add(hists.backgroundAntiOmegaSum[1]);
        for(int j = 3; j < hists.backgroundOmegaSum.size(); j+=2) {
            TH1F* temp = (TH1F*)hists.backgroundOmegaSum[j]->Clone(Form("mergedBkgSumOmegaMixed_%d", j));
            temp->Add(hists.backgroundAntiOmegaSum[j]);
            mergedMixedEventBkgSumVec.push_back(temp);
        }
    for(size_t i = 0; i < hists.backgroundOmegaSum.size()/2; i++){
        TH1F* cf = createCFHistogram(mergedSameEventBkgSumVec[i], mergedMixedEventBkgSumVec[i], Form("cf_merged_bkg_sum_pair_antipair_%zu", i), rebinFactor);
        cfHists.backgroundOmegaSum.push_back(cf);
    }

    return cfHists;
}

std::vector<TH1F*> bkgAvg(HistogramGroups& hists){ // for merged particle antiparticle cf
    std::vector<TH1F*> avgBgHists;

    for(size_t i = 0; i < hists.backgroundLowOmega.size(); i++){
        TH1F* avg = (TH1F*)hists.backgroundLowOmega[i]->Clone(Form("avgBkg_%zu", i));
        avg->Add(hists.backgroundHighOmega[i]);
        avg->Scale(0.5);
        avgBgHists.push_back(avg);
    }
    return avgBgHists;

}


std::vector<TH1F*> cfBackgroundCorrection(TString filename, Int_t rebinFactor=1, Double_t signalBgRatio = 2.08784){
    
    std::vector<TH1F*> correctedCFs;
    HistogramGroups mergedHists = loadkStarHistograms(filename);

    std::cout<< mergedHists.omegaCor.size() << std::endl; // assuming all have the same number of histograms
    HistogramGroups mergedkStarMergedPairAntipair = kStarMergedPairAntipair(mergedHists);   //zmienia histogramy mergedHits

    for(size_t i = 0; i < mergedkStarMergedPairAntipair.backgroundOmegaSum.size(); i++)
        mergedkStarMergedPairAntipair.backgroundOmegaSum[i]->Scale(0.5*signalBgRatio);
    
    std::vector<TH1F*> correctedSameEventVec;
        for(int j = 0; j < mergedkStarMergedPairAntipair.omegaCor.size(); j+=2) {
            TH1F* temp = (TH1F*)mergedkStarMergedPairAntipair.omegaCor[j]->Clone(Form("correctedOmegakStar_%d", j));
            temp->Add(mergedkStarMergedPairAntipair.backgroundOmegaSum[j], -1.0);
            correctedSameEventVec.push_back(temp);
        }

    TCanvas *tc = new TCanvas("Corrected same event", "Corrected same event", 800, 600);
    tc->Divide(2,3);
    for(size_t i = 0; i < correctedSameEventVec.size(); i++)    {
        configureHistogram(correctedSameEventVec[i], kRed, kFullCircle, 0.5);
        tc->cd(i+1);
        correctedSameEventVec[i]->GetXaxis()->SetRangeUser(0,0.5);
        correctedSameEventVec[i]->Draw();
    }
    
    std::vector<TH1F*> correctedMixedEventVec;
    
        for(int j = 1; j < mergedkStarMergedPairAntipair.omegaCor.size(); j+=2) {
            TH1F* temp = (TH1F*)mergedkStarMergedPairAntipair.omegaCor[j]->Clone(Form("correctedOmegakStarMixed_%d", j));
            temp->Add(mergedkStarMergedPairAntipair.backgroundOmegaSum[j], -1.0);
            temp->SetName(Form("correctedOmegakStarMixed_%d", j));
            correctedMixedEventVec.push_back(temp);
        }
    
    for(size_t i = 0; i < mergedkStarMergedPairAntipair.omegaCor.size()/2; i++){
        TH1F* cf = createCFHistogram(correctedSameEventVec[i], correctedMixedEventVec[i], Form("cf_corrected_pair_antipair_%zu", i), rebinFactor);
        correctedCFs.push_back(cf);
    }

    return correctedCFs;
}

void testCorrection(){
    std::vector<TH1F*> correctedCFs = cfBackgroundCorrection("kStarMerged.root",5);
    std::cout << "Number of corrected CFs: " << correctedCFs.size() << std::endl;
    TCanvas* cm = new TCanvas("cor CFs", "Merged CFs", 800, 600);
    cm->Divide(2,3);
    for(size_t i = 0; i < correctedCFs.size(); i++) {
        configureHistogram(correctedCFs[i], kBlue, kFullCircle, 0.5);
        cm->cd(i+1);
        correctedCFs[i]->GetXaxis()->SetRangeUser(0,0.5);
        correctedCFs[i]->GetYaxis()->SetRangeUser(0.5,1.5);
        correctedCFs[i]->Draw();
    }
}

void drawCFs(){


    std::cout<< "Loading histograms from file..." << std::endl;
    HistogramGroups mergedHists = loadkStarHistograms("kStarMergedOld.root");

    std::cout<< mergedHists.omegaCor.size() << std::endl; // assuming all have the same number of histograms

    int rebinFactor = 5;

    // HistogramGroups mergedCFs = kStarToCFs(mergedHists, rebinFactor);                //zmienia histogramy mergedHits
    HistogramGroups mergedCFsMergedPairAntipair = kStarToCFsMergedpairAntipair(mergedHists, rebinFactor);   //zmienia histogramy mergedHits

    std::cout<< "Merged histograms created" << std::endl;
    std::cout<< "number of bins in kstarCF: " << mergedCFsMergedPairAntipair.omegaCor[0]->GetNbinsX() << std::endl;    

    TCanvas* cm = new TCanvas("Merged CFs", "Merged CFs", 800, 600);
    cm->Divide(2,3);
    for(size_t i = 0; i < mergedCFsMergedPairAntipair.omegaCor.size(); i++) {
        configureHistogram(mergedCFsMergedPairAntipair.omegaCor[i], kBlue, kFullCircle, 0.5);
        cm->cd(i+1);
        mergedCFsMergedPairAntipair.omegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFsMergedPairAntipair.omegaCor[i]->Draw();
    }

    TCanvas* cm2 = new TCanvas("Merged CFs background", "Merged CFs background", 800, 600);
    cm2->Divide(2,3);
    for(size_t i = 0; i < mergedCFsMergedPairAntipair.backgroundOmegaSum.size(); i++) {
        configureHistogram(mergedCFsMergedPairAntipair.backgroundOmegaSum[i], kBlue, kFullCircle, 0.5);
        cm2->cd(i+1);
        mergedCFsMergedPairAntipair.backgroundOmegaSum[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFsMergedPairAntipair.backgroundOmegaSum[i]->Draw();
    }

    TCanvas* cm3 = new TCanvas("Merged CFs background high", "Merged CFs background high", 800, 600);
    cm3->Divide(2,3);
    for(size_t i = 0; i < mergedCFsMergedPairAntipair.backgroundHighOmega.size(); i++) {
        configureHistogram(mergedCFsMergedPairAntipair.backgroundHighOmega[i], kBlue, kFullCircle, 0.5);
        cm3->cd(i+1);
        mergedCFsMergedPairAntipair.backgroundHighOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFsMergedPairAntipair.backgroundHighOmega[i]->Draw();
    }

    TCanvas* cm4 = new TCanvas("Merged CFs background low", "Merged CFs background low", 800, 600);
    cm4->Divide(2,3);
    for(size_t i = 0; i < mergedCFsMergedPairAntipair.backgroundLowOmega.size(); i++) {
        configureHistogram(mergedCFsMergedPairAntipair.backgroundLowOmega[i], kBlue, kFullCircle, 0.5);
        cm4->cd(i+1);
        mergedCFsMergedPairAntipair.backgroundLowOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFsMergedPairAntipair.backgroundLowOmega[i]->Draw();
    }

    TCanvas* cm5 = new TCanvas("Merged CFs signal + bkg", "Merged CFs signal + bkg", 800, 600);
    cm5->Divide(2,3);
    for(size_t i = 0; i < mergedCFsMergedPairAntipair.omegaCor.size(); i++) {
        configureHistogram(mergedCFsMergedPairAntipair.omegaCor[i], kBlue, kFullCircle, 0.75);
        configureHistogram(mergedCFsMergedPairAntipair.backgroundHighOmega[i], kRed, kFullSquare, 0.75);
        configureHistogram(mergedCFsMergedPairAntipair.backgroundLowOmega[i], kGreen, kFullDiamond, 0.75);
        cm5->cd(i+1);
        mergedCFsMergedPairAntipair.omegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
        mergedCFsMergedPairAntipair.omegaCor[i]->Draw();
        mergedCFsMergedPairAntipair.backgroundHighOmega[i]->Draw("SAME");
        mergedCFsMergedPairAntipair.backgroundLowOmega[i]->Draw("SAME");
    }
    TLegend* legend = new TLegend();
    legend->AddEntry(mergedCFsMergedPairAntipair.omegaCor[0],"Merged CF","p");
    legend->AddEntry(mergedCFsMergedPairAntipair.backgroundHighOmega[0],"Merged high bkg CF","p");
    legend->AddEntry(mergedCFsMergedPairAntipair.backgroundLowOmega[0],"Merged low bkg CF","p");
    cm5->cd(1);
    legend->Draw();

    TCanvas* cm6 = new TCanvas("Cf comparison", "CF comparison", 800, 600);
    cm6->cd();
    mergedCFsMergedPairAntipair.omegaCor[1]->Draw();
    mergedCFsMergedPairAntipair.omegaCor[1]->GetYaxis()->SetTitle("C(k*)");
    mergedCFsMergedPairAntipair.omegaCor[1]->SetStats(0);
    mergedCFsMergedPairAntipair.omegaCor[2]->SetStats(0);
    mergedCFsMergedPairAntipair.omegaCor[3]->SetStats(0);
    mergedCFsMergedPairAntipair.omegaCor[4]->SetStats(0);
    mergedCFsMergedPairAntipair.omegaCor[5]->SetStats(0);
    configureHistogram(mergedCFsMergedPairAntipair.omegaCor[1], kBlue, kFullCircle, 1.2);
    mergedCFsMergedPairAntipair.omegaCor[2]->Draw("SAME");
    configureHistogram(mergedCFsMergedPairAntipair.omegaCor[2], kRed, kFullSquare, 1.2);
    mergedCFsMergedPairAntipair.omegaCor[3]->Draw("SAME");
    configureHistogram(mergedCFsMergedPairAntipair.omegaCor[3], 8, kFullDiamond, 1.2);
    mergedCFsMergedPairAntipair.omegaCor[4]->Draw("SAME");
    configureHistogram(mergedCFsMergedPairAntipair.omegaCor[4], 9, kFullCross, 1.2);
    mergedCFsMergedPairAntipair.omegaCor[5]->Draw("SAME");
    configureHistogram(mergedCFsMergedPairAntipair.omegaCor[5], 1, kFullDoubleDiamond, 1.5);
    TLegend* legend2 = new TLegend();
    legend2->AddEntry(mergedCFsMergedPairAntipair.omegaCor[1],"0-10\% centrality","p");
    legend2->AddEntry(mergedCFsMergedPairAntipair.omegaCor[2],"10-20\% centrality","p");
    legend2->AddEntry(mergedCFsMergedPairAntipair.omegaCor[3],"20-30\% centrality","p");
    legend2->AddEntry(mergedCFsMergedPairAntipair.omegaCor[4],"30-50\% centrality","p");
    legend2->AddEntry(mergedCFsMergedPairAntipair.omegaCor[5],"50-100\% centrality","p");
    cm6->cd();
    legend2->Draw();    
    


    // TCanvas *c1 = new TCanvas("pbarOmega","pbarOmega",800,600);
    // c1->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(mergedCFs.omegaCor[i],kBlue,kFullCircle,0.5);
    //     c1->cd(i+1);
    //     mergedCFs.omegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
    //     mergedCFs.omegaCor[i]->Draw();
    // }

    // TCanvas *c2 = new TCanvas("pOmegabar","pOmegabar",800,600);
    // c2->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(mergedCFs.antiOmegaCor[i],kRed,kFullSquare,0.5);
    //     c2->cd(i+1);
    //     mergedCFs.antiOmegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
    //     mergedCFs.antiOmegaCor[i]->Draw();
    // }    



    // TCanvas *c10 = new TCanvas("Background high","Background high",800,600);
    // c10->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(mergedCFs.backgroundHighOmega[i], kBlue, kFullCircle, 0.5);
    //     configureHistogram(mergedCFs.backgroundHighAntiOmega[i], kRed, kFullSquare, 0.5);
    //     mergedCFs.backgroundHighOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
    //     c10->cd(i+1);
    //     mergedCFs.backgroundHighOmega[i]->Draw();
    //     mergedCFs.backgroundHighAntiOmega[i]->Draw("SAME");
    // }
    // c10->cd(1);
    // TLegend *legendHigh = new TLegend();
    // legendHigh->AddEntry(mergedCFs.backgroundHighOmega[0],"Omega background","p");
    // legendHigh->AddEntry(mergedCFs.backgroundHighAntiOmega[0],"Antyomega background","p");
    // legendHigh->Draw();


    // TCanvas *c11 = new TCanvas("Background low","Background low",800,600);
    // c11->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(mergedCFs.backgroundLowOmega[i], kBlue, kFullCircle, 0.5);
    //     configureHistogram(mergedCFs.backgroundLowAntiOmega[i], kRed, kFullSquare, 0.5);
    //     mergedCFs.backgroundLowOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
    //     c11->cd(i+1);
    //     mergedCFs.backgroundLowOmega[i]->Draw();
    //     mergedCFs.backgroundLowAntiOmega[i]->Draw("SAME");
    // }
    // c11->cd(1);
    // TLegend *legendLow = new TLegend();
    // legendLow->AddEntry(mergedCFs.backgroundLowOmega[0],"Omega background","p");
    // legendLow->AddEntry(mergedCFs.backgroundLowAntiOmega[0],"Antyomega background","p");
    // legendLow->Draw();



    // TCanvas *c5 = new TCanvas("Cf comparison","Cf comparison",800,600);
    // c5->Divide(2,3);
    // for(Int_t i = 0; i < 6; i++){
    //     c5->cd(i+1);
    //     mergedCFs.omegaCor[i]->Draw();
    //     mergedCFs.antiOmegaCor[i]->Draw("SAME");
    // }

    // TLegend *legend3 = new TLegend();
    // legend3->AddEntry(mergedCFs.omegaCor[0],"antiproton-omega","p");
    // legend3->AddEntry(mergedCFs.antiOmegaCor[0],"proton-antiomega","p");
    // c5->cd(1);
    // legend3->Draw();



}