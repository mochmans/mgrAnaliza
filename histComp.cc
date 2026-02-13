
#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <string>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include "TLegend.h"
#include "TDirectory.h"
#include "cfHistRead.cc"

void configureHistogram(TH1* histo, Color_t markerColor, Style_t markerStyle, Size_t markerSize){
    histo->SetMarkerColor(markerColor);
    histo->SetMarkerStyle(markerStyle);
    histo->SetMarkerSize(markerSize);
    histo->GetXaxis()->SetRangeUser(0,1.5);
    histo->Scale(1.0);
}


int histComp(){

    TH1::SetDefaultSumw2();
    
    HistogramGroups zzhHists = loadCorrelationHistograms("zzh_medium_correlationResultsFull.root");
    HistogramGroups zzoHists = loadCorrelationHistograms("zzo_medium_correlationResultsFull.root");
    HistogramGroups zzmHists = loadCorrelationHistograms("zzm_medium_correlationResultsFull.root");

    double_t markerSize = 0.7;

    if(!(zzhHists.omegaCor.empty() || zzoHists.omegaCor.empty() || zzmHists.omegaCor.empty())) {
        TCanvas *c1 = new TCanvas("Correlation functions comparison omega","Correlation functions comparison omega",800,600);
        c1->Divide(2,3);
        for(Int_t i = 0; i<6; i++){
            configureHistogram(zzhHists.omegaCor[i], kBlue, kFullCircle, markerSize );
            configureHistogram(zzoHists.omegaCor[i], kRed, kFullSquare, markerSize );
            configureHistogram(zzmHists.omegaCor[i], kGreen, kFullDiamond, markerSize );
            c1->cd(i+1);
            zzhHists.omegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
            zzhHists.omegaCor[i]->GetYaxis()->SetRangeUser(0.0,2.5);
            zzhHists.omegaCor[i]->Draw();
            zzoHists.omegaCor[i]->Draw("SAME");
            zzmHists.omegaCor[i]->Draw("SAME");
        }
        c1->cd(1);
        TLegend *l1 = new TLegend();
            l1->AddEntry(zzhHists.omegaCor[0],"ZZH","p");
            l1->AddEntry(zzoHists.omegaCor[0],"ZZO","p");
            l1->AddEntry(zzmHists.omegaCor[0],"ZZM","p");
            l1->Draw();
    }

    if(!(zzhHists.antiOmegaCor.empty() || zzoHists.antiOmegaCor.empty() || zzmHists.antiOmegaCor.empty())) {
        TCanvas *c2 = new TCanvas("Correlation functions comparison antiomega","Correlation functions comparison antiomega",800,600);
        c2->Divide(2,3);
        for(Int_t i = 0; i<6; i++){
            configureHistogram(zzhHists.antiOmegaCor[i], kBlue, kFullCircle, markerSize );
            configureHistogram(zzoHists.antiOmegaCor[i], kRed, kFullSquare, markerSize );
            configureHistogram(zzmHists.antiOmegaCor[i], kGreen, kFullTriangleUp, markerSize );
            c2->cd(i+1);
            zzhHists.antiOmegaCor[i]->GetXaxis()->SetRangeUser(0,0.5);
            zzhHists.antiOmegaCor[i]->GetYaxis()->SetRangeUser(0.0,2.5);
            zzhHists.antiOmegaCor[i]->Draw();
            zzoHists.antiOmegaCor[i]->Draw("SAME");
            zzmHists.antiOmegaCor[i]->Draw("SAME");
        }
        c2->cd(1);
        TLegend *l2 = new TLegend();
            l2->AddEntry(zzhHists.antiOmegaCor[0],"ZZH","p");
            l2->AddEntry(zzoHists.antiOmegaCor[0],"ZZO","p");
            l2->AddEntry(zzmHists.antiOmegaCor[0],"ZZM","p");
            l2->Draw();
    }
    
    if(!(zzhHists.backgroundOmegaSum.empty() || zzoHists.backgroundOmegaSum.empty() || zzmHists.backgroundOmegaSum.empty())) {
        TCanvas *c3 = new TCanvas("Correlation functions comparison background sum","Correlation functions comparison background sum",800,600);
        c3->Divide(2,3);
        for(Int_t i = 0; i<6; i++){
            configureHistogram(zzhHists.backgroundOmegaSum[i], kBlue, kFullCircle, markerSize );
            configureHistogram(zzoHists.backgroundOmegaSum[i], kRed, kFullSquare, markerSize );
            configureHistogram(zzmHists.backgroundOmegaSum[i], kGreen, kFullTriangleUp, markerSize );
            c3->cd(i+1);
            zzhHists.backgroundOmegaSum[i]->GetXaxis()->SetRangeUser(0,0.5);
            zzhHists.backgroundOmegaSum[i]->GetYaxis()->SetRangeUser(0.0,2.5);
            zzhHists.backgroundOmegaSum[i]->Draw();
            zzoHists.backgroundOmegaSum[i]->Draw("SAME");
            zzmHists.backgroundOmegaSum[i]->Draw("SAME");
        }
        c3->cd(1);
        TLegend *l3 = new TLegend();
            l3->AddEntry(zzhHists.backgroundOmegaSum[0],"ZZH","p");
            l3->AddEntry(zzoHists.backgroundOmegaSum[0],"ZZO","p");
            l3->AddEntry(zzmHists.backgroundOmegaSum[0],"ZZM","p");
            l3->Draw(); 
    }

    if(!(zzhHists.backgroundHighOmega.empty() || zzoHists.backgroundHighOmega.empty() || zzmHists.backgroundHighOmega.empty())) {
    TCanvas *c4 = new TCanvas("Correlation functions comparison background high","Correlation functions comparison background high",800,600);
    c4->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(zzhHists.backgroundHighOmega[i], kBlue, kFullCircle, markerSize );
        configureHistogram(zzoHists.backgroundHighOmega[i], kRed, kFullSquare, markerSize );
        configureHistogram(zzmHists.backgroundHighOmega[i], kGreen, kFullTriangleUp, markerSize );
        c4->cd(i+1);
        zzhHists.backgroundHighOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        zzhHists.backgroundHighOmega[i]->GetYaxis()->SetRangeUser(0.0,2.5);
        zzhHists.backgroundHighOmega[i]->Draw();
        zzoHists.backgroundHighOmega[i]->Draw("SAME");
        zzmHists.backgroundHighOmega[i]->Draw("SAME");
    }
    c4->cd(1);
    TLegend *l4 = new TLegend();
        l4->AddEntry(zzhHists.backgroundHighOmega[0],"ZZH","p");
        l4->AddEntry(zzoHists.backgroundHighOmega[0],"ZZO","p");
        l4->AddEntry(zzmHists.backgroundHighOmega[0],"ZZM","p");
        l4->Draw(); 
    }
    
    if(!(zzhHists.backgroundLowOmega.empty() || zzoHists.backgroundLowOmega.empty() || zzmHists.backgroundLowOmega.empty())) {  
    TCanvas *c5 = new TCanvas("Correlation functions comparison background low","Correlation functions comparison background low",800,600);
    c5->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(zzhHists.backgroundLowOmega[i], kBlue, kFullCircle   , markerSize );
        configureHistogram(zzoHists.backgroundLowOmega[i], kRed, kFullSquare, markerSize );
        configureHistogram(zzmHists.backgroundLowOmega[i], kGreen, kFullTriangleUp, markerSize );
        c5->cd(i+1);
        zzhHists.backgroundLowOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        zzhHists.backgroundLowOmega[i]->GetYaxis()->SetRangeUser(0.0,2.5);
        zzhHists.backgroundLowOmega[i]->Draw();
        zzoHists.backgroundLowOmega[i]->Draw("SAME");
        zzmHists.backgroundLowOmega[i]->Draw("SAME");
    }
    c5->cd(1);
    TLegend *l5 = new TLegend();
        l5->AddEntry(zzhHists.backgroundLowOmega[0],"ZZH","p");
        l5->AddEntry(zzoHists.backgroundLowOmega[0],"ZZO","p");
        l5->AddEntry(zzmHists.backgroundLowOmega[0],"ZZM","p");
        l5->Draw(); 
    }      
    
    if(!(zzhHists.backgroundAntiOmegaSum.empty() || zzoHists.backgroundAntiOmegaSum.empty() || zzmHists.backgroundAntiOmegaSum.empty())) {
    TCanvas *c6 = new TCanvas("Correlation functions comparison background antiomega sum","Correlation functions comparison background antiomega sum",800,600); 
    c6->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(zzhHists.backgroundAntiOmegaSum[i], kBlue, kFullCircle, markerSize );
        configureHistogram(zzoHists.backgroundAntiOmegaSum[i], kRed, kFullSquare, markerSize );
        configureHistogram(zzmHists.backgroundAntiOmegaSum[i], kGreen, kFullTriangleUp, markerSize );
        c6->cd(i+1);
        zzhHists.backgroundAntiOmegaSum[i]->GetXaxis()->SetRangeUser(0,0.5);
        zzhHists.backgroundAntiOmegaSum[i]->GetYaxis()->SetRangeUser(0.0,2.5);
        zzhHists.backgroundAntiOmegaSum[i]->Draw();
        zzoHists.backgroundAntiOmegaSum[i]->Draw("SAME");
        zzmHists.backgroundAntiOmegaSum[i]->Draw("SAME");   
    }
    c6->cd(1);
    TLegend *l6 = new TLegend();    
        l6->AddEntry(zzhHists.backgroundAntiOmegaSum[0],"ZZH","p");
        l6->AddEntry(zzoHists.backgroundAntiOmegaSum[0],"ZZO","p");
        l6->AddEntry(zzmHists.backgroundAntiOmegaSum[0],"ZZM","p");
        l6->Draw(); 
    }

    if(!(zzhHists.backgroundHighAntiOmega.empty() || zzoHists.backgroundHighAntiOmega.empty() || zzmHists.backgroundHighAntiOmega.empty())) {
    TCanvas *c7 = new TCanvas("Correlation functions comparison background antiomega high","Correlation functions comparison background antiomega high",800,600);
    c7->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(zzhHists.backgroundHighAntiOmega[i], kBlue, kFullCircle, markerSize );
        configureHistogram(zzoHists.backgroundHighAntiOmega[i], kRed, kFullSquare, markerSize );
        configureHistogram(zzmHists.backgroundHighAntiOmega[i], kGreen, kFullTriangleUp, markerSize );
        c7->cd(i+1);
        zzhHists.backgroundHighAntiOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        zzhHists.backgroundHighAntiOmega[i]->GetYaxis()->SetRangeUser(0.0,2.5);
        zzhHists.backgroundHighAntiOmega[i]->Draw();
        zzoHists.backgroundHighAntiOmega[i]->Draw("SAME");
        zzmHists.backgroundHighAntiOmega[i]->Draw("SAME");
    }
    c7->cd(1);
    TLegend *l7 = new TLegend();
        l7->AddEntry(zzhHists.backgroundHighAntiOmega[0],"ZZH","p");
        l7->AddEntry(zzoHists.backgroundHighAntiOmega[0],"ZZO","p");
        l7->AddEntry(zzmHists.backgroundHighAntiOmega[0],"ZZM","p");
        l7->Draw(); 
    }
    
    if(!(zzhHists.backgroundLowAntiOmega.empty() || zzoHists.backgroundLowAntiOmega.empty() || zzmHists.backgroundLowAntiOmega.empty())) {
    TCanvas *c8 = new TCanvas("Correlation functions comparison background antiomega low","Correlation functions comparison background antiomega low",800,600);
    c8->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(zzhHists.backgroundLowAntiOmega[i], kBlue, kFullCircle, markerSize );
        configureHistogram(zzoHists.backgroundLowAntiOmega[i], kRed, kFullSquare, markerSize );
        configureHistogram(zzmHists.backgroundLowAntiOmega[i], kGreen, kFullTriangleUp, markerSize );
        c8->cd(i+1);
        zzhHists.backgroundLowAntiOmega[i]->GetXaxis()->SetRangeUser(0,0.5);
        zzhHists.backgroundLowAntiOmega[i]->GetYaxis()->SetRangeUser(0.0,2.5);
        zzhHists.backgroundLowAntiOmega[i]->Draw();
        zzoHists.backgroundLowAntiOmega[i]->Draw("SAME");
        zzmHists.backgroundLowAntiOmega[i]->Draw("SAME");
    }
    c8->cd(1);
    TLegend *l8 = new TLegend();
        l8->AddEntry(zzhHists.backgroundLowAntiOmega[0],"ZZH","p");
        l8->AddEntry(zzoHists.backgroundLowAntiOmega[0],"ZZO","p");
        l8->AddEntry(zzmHists.backgroundLowAntiOmega[0],"ZZM","p");
        l8->Draw(); 
    }
    
    return 0;
}

