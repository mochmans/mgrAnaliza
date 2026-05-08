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

std::vector<TF1*> backgroundFit(TString filename, TString directoryName, bool drawFits = kFALSE) {

    TFile* f = new TFile(filename,"open");
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir = (TDirectory*)f->Get(directoryName);

    TDirectory *Cascade = (TDirectory*)Dir->Get("Cascade");

    TH1D *hInvMassOmega = (TH1D*)Cascade->Get("hInvMassOmega");

    std::cout << "minv integral in 1.66358,1.666965 range: " << hInvMassOmega->Integral(hInvMassOmega->GetXaxis()->FindBin(1.66358), hInvMassOmega->GetXaxis()->FindBin(1.666965)) << std::endl;

    TF1 *fitFunc = new TF1("fitFunc","gaus(0) + gaus(3) + pol5(6)",1.65,1.70);
    // double pars[12]={114269,1.67163,0.00157569,29456.1,1.67161,0.00359362,-6.4086e6,739425,-1.83344e6,1.10416e6,170405,-487243};
    double pars[12] = {120000,1.6715,0.001,60000,1.6715,0.0025,80000};
    fitFunc->SetParameters(pars);
    // fitFunc->SetParLimits(2,0.0003,0.003);
    // fitFunc->SetParLimits(5,0.001,0.01);
    hInvMassOmega->Fit(fitFunc,"R");   
    hInvMassOmega->GetXaxis()->SetRangeUser(1.65,1.70);
    // hInvMassOmega->GetYaxis()->SetRangeUser(4500,8500);
    hInvMassOmega->Draw();
  
    TF1* gaussFit1 = new TF1("gaussFit1", "gaus", 1.65, 1.70);
for(int i=0; i<gaussFit1->GetNpar(); i++) 
    gaussFit1->SetParameter(i, fitFunc->GetParameter(i));

    TF1* gaussFit2 = new TF1("gaussFit2", "gaus", 1.65, 1.70);
for(int i=0; i<gaussFit2->GetNpar(); i++) 
    gaussFit2->SetParameter(i, fitFunc->GetParameter(i + gaussFit1->GetNpar()));    



TF1* pol4Fit = new TF1("pol4Fit", "pol5", 1.65, 1.70);
for(int i=0; i<pol4Fit->GetNpar() ; i++) 
    pol4Fit->SetParameter(i, fitFunc->GetParameter(i + gaussFit1->GetNpar() + gaussFit2->GetNpar()));

double_t integralLowBg = pol4Fit->Integral(1.66358,1.666965);
double_t integralSignal = pol4Fit->Integral(1.66865,1.67525);
double_t integralHighBg = pol4Fit->Integral(1.677,1.68);
double_t integralSignalPeak = gaussFit1->Integral(1.66865,1.67525);
integralSignalPeak += gaussFit2->Integral(1.66865,1.67525);

std::cout << "pol2 integral in range 1.66358 to 1.666965: " << integralLowBg << " in dx of " << 1.666965 - 1.66358 << std::endl;
std::cout << "pol2 integral in range 1.66865 to 1.67525 (signal): "<< integralSignal << " in dx of " << 1.67525 - 1.66865 << std::endl;
std::cout << "pol2 integral in range 1.677 to 1.68: " << integralHighBg << " in dx of " << 1.68 - 1.677 << std::endl;
std::cout << "gauss integral in range 1.66865 to 1.67525 (signal peak): " << integralSignalPeak << " in dx of " << 1.67525 - 1.66865 << std::endl;
std::cout<< "Signal purity =" << integralSignalPeak / (integralSignalPeak + integralSignal) << std::endl;

double_t integralAvgBg = (integralLowBg + integralHighBg) / 2.0;
std::cout << "Average background integral: " << integralAvgBg << std::endl;
std::cout << "Signal to background ratio (using average background): " << integralSignal / integralAvgBg << std::endl;

if(drawFits) {
    TCanvas *c1 = new TCanvas("Invariant Mass Omega","Invariant Mass Omega",800,600);
    hInvMassOmega->Draw();
    fitFunc->Draw("SAME");
    gaussFit1->SetLineColor(kBlue+4);
    gaussFit1->Draw("SAME");
    gaussFit2->SetLineColor(kGreen+3);
    gaussFit2->Draw("SAME");
    pol4Fit->SetLineColor(kGreen);
    pol4Fit->Draw("SAME");
}

    std::vector<TF1*> fitFunctions;
    fitFunctions.push_back(fitFunc);
    fitFunctions.push_back(gaussFit1);
    fitFunctions.push_back(gaussFit2);
    fitFunctions.push_back(pol4Fit);
    return fitFunctions;
}