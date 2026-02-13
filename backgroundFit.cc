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

std::vector<TF1*> backgroundFit(TFile* filename, TString directoryName, bool drawFits = kFALSE) {

    TFile* f = filename;
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir = (TDirectory*)f->Get(directoryName);

    TDirectory *Cascade = (TDirectory*)Dir->Get("Cascade");

    TH1D *hInvMassOmega = (TH1D*)Cascade->Get("hInvMassOmega");

    TF1 *fitFunc = new TF1("fitFunc","gaus(0)+pol3(3)",1.65,1.70);
    fitFunc->SetParameters(1,1.67,0.02,1,1,1);
    hInvMassOmega->Fit(fitFunc,"R");   
    hInvMassOmega->GetXaxis()->SetRangeUser(1.65,1.70);
    // hInvMassOmega->GetYaxis()->SetRangeUser(4500,8500);
    hInvMassOmega->Draw();
  
    TF1* gaussFit = new TF1("gaussFit", "gaus", 1.65, 1.70);
for(int i=0; i<gaussFit->GetNpar(); i++) 
    gaussFit->SetParameter(i, fitFunc->GetParameter(i));


TF1* pol3Fit = new TF1("pol3Fit", "pol3", 1.65, 1.70);
for(int i=0; i<pol3Fit->GetNpar() ; i++) 
    pol3Fit->SetParameter(i, fitFunc->GetParameter(i + gaussFit->GetNpar()));

if(drawFits) {
    TCanvas *c1 = new TCanvas("Invariant Mass Omega","Invariant Mass Omega",800,600);
    hInvMassOmega->Draw();
    fitFunc->Draw("SAME");
    gaussFit->SetLineColor(kBlue);
    gaussFit->Draw("SAME");
    pol3Fit->SetLineColor(kGreen);
    pol3Fit->Draw("SAME");
}

    std::vector<TF1*> fitFunctions;
    fitFunctions.push_back(fitFunc);
    fitFunctions.push_back(gaussFit);
    fitFunctions.push_back(pol3Fit);
    return fitFunctions;
}