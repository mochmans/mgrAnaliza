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



void configureHistogram(TH1* histo, Color_t markerColor, Style_t markerStyle, Size_t markerSize){
    histo->SetMarkerColor(markerColor);
    histo->SetMarkerStyle(markerStyle);
    histo->SetMarkerSize(markerSize);
    histo->GetXaxis()->SetRangeUser(0,1.5);
    histo->Scale(1.0);
}



int ratioFull(){
    TH1::SetDefaultSumw2();
    // TString rootFileName = "mediumDataOmegaARfull.root";
    TString rootFileName = "zzo_medium.root";
    TFile* f = new TFile(rootFileName,"open");
    // TFile* fFit = new TFile("mediumDataOmegaARWide.root","open");

    Int_t rebinFactor = 1;
    
    std::vector<TH1F*> backgroundHighOmega = correlationFunctions(f, "femto-universe-pair-task-track-cascade-extended_omega_bkg_high", rebinFactor);
    std::vector<TH1F*> backgroundLowOmega = correlationFunctions(f, "femto-universe-pair-task-track-cascade-extended_omega_bkg_low", rebinFactor);
    std::vector<TH1F*> backgroundHighAntiOmega = correlationFunctions(f, "femto-universe-pair-task-track-cascade-extended_antyomega_bkg_high", rebinFactor);
    std::vector<TH1F*> backgroundLowAntiOmega = correlationFunctions(f, "femto-universe-pair-task-track-cascade-extended_antyomega_bkg_low", rebinFactor);

    std::vector<TH1F*> backgroundOmegaSum = backgroundCorrelationSum(f, "femto-universe-pair-task-track-cascade-extended_omega_bkg_low", "femto-universe-pair-task-track-cascade-extended_omega_bkg_high", rebinFactor);
    std::vector<TH1F*> backgroundAntiOmegaSum = backgroundCorrelationSum(f, "femto-universe-pair-task-track-cascade-extended_antyomega_bkg_low", "femto-universe-pair-task-track-cascade-extended_antyomega_bkg_high", rebinFactor);

    std::vector<TH1F*> omegaCorFuncs = correlationFunctions(f,"femto-universe-pair-task-track-cascade-extended",rebinFactor);
    std::vector<TH1F*> antiOmegaCorFuncs = correlationFunctions(f,"femto-universe-pair-task-track-cascade-extended_antyomega",rebinFactor);

    

    
    
    TCanvas *c12 = new TCanvas("Background sum low+high","Background sum low+high",800,600);
    c12->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(backgroundOmegaSum[i], kBlue, kFullCircle, 0.5);
        configureHistogram(backgroundAntiOmegaSum[i], kRed, kFullSquare, 0.5);
        backgroundOmegaSum[i]->GetXaxis()->SetRangeUser(0,0.5);
        c12->cd(i+1);
        backgroundOmegaSum[i]->Draw();
        backgroundAntiOmegaSum[i]->Draw("SAME");
    }
    TLegend *legendSum = new TLegend();
    legendSum->AddEntry(backgroundOmegaSum[0],"Omega background","p");
    legendSum->AddEntry(backgroundAntiOmegaSum[0],"Antyomega background","p");
    c12->cd(1);
    legendSum->Draw();

    TCanvas *c1 = new TCanvas("pbarOmega","pbarOmega",800,600);
    c1->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(omegaCorFuncs[i],kBlue,kFullCircle,0.5);
        c1->cd(i+1);
        omegaCorFuncs[i]->GetXaxis()->SetRangeUser(0,0.5);
        omegaCorFuncs[i]->Draw();
    }


    TCanvas *c2 = new TCanvas("pOmegabar","pOmegabar",800,600);
    c2->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(antiOmegaCorFuncs[i],kRed,kFullSquare,0.5);
        c2->cd(i+1);
        antiOmegaCorFuncs[i]->GetXaxis()->SetRangeUser(0,0.5);
        antiOmegaCorFuncs[i]->Draw();
    }    



    TCanvas *c10 = new TCanvas("Background high","Background high",800,600);
    c10->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(backgroundHighOmega[i], kBlue, kFullCircle, 0.5);
        configureHistogram(backgroundHighAntiOmega[i], kRed, kFullSquare, 0.5);
        backgroundHighOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        c10->cd(i+1);
        backgroundHighOmega[i]->Draw();
        backgroundHighAntiOmega[i]->Draw("SAME");
    }
    c10->cd(1);
    TLegend *legendHigh = new TLegend();
    legendHigh->AddEntry(backgroundHighOmega[0],"Omega background","p");
    legendHigh->AddEntry(backgroundHighAntiOmega[0],"Antyomega background","p");
    legendHigh->Draw();


    TCanvas *c11 = new TCanvas("Background low","Background low",800,600);
    c11->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(backgroundLowOmega[i], kBlue, kFullCircle, 0.5);
        configureHistogram(backgroundLowAntiOmega[i], kRed, kFullSquare, 0.5);
        backgroundLowOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        c11->cd(i+1);
        backgroundLowOmega[i]->Draw();
        backgroundLowAntiOmega[i]->Draw("SAME");
    }
    c11->cd(1);
    TLegend *legendLow = new TLegend();
    legendLow->AddEntry(backgroundLowOmega[0],"Omega background","p");
    legendLow->AddEntry(backgroundLowAntiOmega[0],"Antyomega background","p");
    legendLow->Draw();



    TCanvas *c5 = new TCanvas("Cf comparison","Cf comparison",800,600);
    c5->Divide(2,3);
    for(Int_t i = 0; i < 6; i++){
        c5->cd(i+1);
        omegaCorFuncs[i]->Draw();
        antiOmegaCorFuncs[i]->Draw("SAME");
    }

    TLegend *legend3 = new TLegend();
    legend3->AddEntry(omegaCorFuncs[0],"antiproton-omega","p");
    legend3->AddEntry(antiOmegaCorFuncs[0],"proton-antiomega","p");
    c5->cd(1);
    legend3->Draw();

    TString resultsFileName = rootFileName.ReplaceAll(".root","_correlationResultsFull.root");

    SaveCorrelationHistograms(
        resultsFileName,
        omegaCorFuncs,
        antiOmegaCorFuncs
        // ,backgroundOmegaSum,
        // backgroundAntiOmegaSum,
        // backgroundHighOmega,
        // backgroundHighAntiOmega,
        // backgroundLowOmega,
        // backgroundLowAntiOmega
    );
     
    return 0; 
}