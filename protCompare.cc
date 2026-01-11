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


int protCompare(){
    TH1::SetDefaultSumw2();

    TFile* f = new TFile("mediumDataOmegaAR2.root","open");

    TDirectory *pbarOmega = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended"); //antyproton omega
    TDirectory *pOmegaBar = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_antyomega"); //proton antyomega
    TDirectory *protons = (TDirectory*)pOmegaBar->Get("Tracks_pos");
    TDirectory *antyprotons = (TDirectory*)pbarOmega->Get("Tracks_neg");

    TH2F *hProtonPt = (TH2F*)protons->Get("hPt");
    TH2F *hAntyProtonPt = (TH2F*)antyprotons->Get("hPt");   

    TH2F *hProtonEta = (TH2F*)protons->Get("hEta");
    TH2F *hAntyProtonEta = (TH2F*)antyprotons->Get("hEta");     
    TH2F *hProtonPhi = (TH2F*)protons->Get("hPhi");
    TH2F *hAntyProtonPhi = (TH2F*)antyprotons->Get("hPhi");


    TH2F *hProtonnSigmaTpc = (TH2F*)protons->Get("nSigmaTPC");
    TH2F *hAntyProtonnSigmaTpc = (TH2F*)antyprotons->Get("nSigmaTPC");
    TH2F *hProtonnSigmaTof = (TH2F*)protons->Get("nSigmaTOF");
    TH2F *hAntyProtonnSigmaTof = (TH2F*)antyprotons->Get("nSigmaTOF");



    TCanvas *c1 = new TCanvas("Protons vs AntyProtons","pvpbar",800,600);
    c1->Divide(3,1);
    c1->cd(1);
    hProtonPt->SetLineColor(kBlue);
    hAntyProtonPt->SetLineColor(kRed);
    hProtonPt->Draw();
    hAntyProtonPt->Draw("SAME");
    c1->cd(2);
    hProtonEta->SetLineColor(kBlue);
    hAntyProtonEta->SetLineColor(kRed);
    hProtonEta->Draw();
    hAntyProtonEta->Draw("SAME");
    c1->cd(3);
    hProtonPhi->SetLineColor(kBlue);
    hAntyProtonPhi->SetLineColor(kRed);
    hProtonPhi->Draw();
    hAntyProtonPhi->Draw("SAME");

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(hProtonPt,"Protons","l");
    leg->AddEntry(hAntyProtonPt,"AntyProtons","l");
    c1->cd(1);
    leg->Draw();

    TCanvas *c2 = new TCanvas("Protons vs AntyProtons nSigma","pvpbar nSigma",800,600);
    c2->Divide(2,2);
    c2->cd(1);
    hProtonnSigmaTpc->Draw();
    c2->cd(2);
    hProtonnSigmaTof->Draw();
    c2->cd(3);
    hAntyProtonnSigmaTpc->Draw();
    c2->cd(4);
    hAntyProtonnSigmaTof->Draw(); 

    return 0;
}