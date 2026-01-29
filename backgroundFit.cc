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

void backgroundFit(){

    TFile* f = new TFile("mediumDataOmegaARWide.root", "READ");
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_omega_wide");

    TDirectory *Cascade = (TDirectory*)Dir->Get("Cascade");

    TH1D *hInvMassOmega = (TH1D*)Cascade->Get("hInvMassOmega");

    TF1 *fitFunc = new TF1("fitFunc","gaus(0)+pol4(3)",1.63,1.71);
    fitFunc->SetParameters(100,1.67,0.02,100,1,0);
    hInvMassOmega->Fit(fitFunc,"R");   
    hInvMassOmega->GetXaxis()->SetRangeUser(1.63,1.71);
    // hInvMassOmega->GetYaxis()->SetRangeUser(4500,8500);
    hInvMassOmega->Draw();
  
    TF1* gaussFit = new TF1("gaussFit", "gaus", 1.63, 1.71);
for(int i=0; i<gaussFit->GetNpar(); i++) 
    gaussFit->SetParameter(i, fitFunc->GetParameter(i));


TF1* pol4Fit = new TF1("pol4Fit", "pol4", 1.63, 1.71);
for(int i=0; i<pol4Fit->GetNpar() ; i++) 
    pol4Fit->SetParameter(i, fitFunc->GetParameter(i + gaussFit->GetNpar()));


    TCanvas *c1 = new TCanvas("Invariant Mass Omega","Invariant Mass Omega",800,600);
    hInvMassOmega->Draw();
    fitFunc->Draw("SAME");
    gaussFit->SetLineColor(kBlue);
    gaussFit->Draw("SAME");
    pol4Fit->SetLineColor(kGreen);
    pol4Fit->Draw("SAME");
}