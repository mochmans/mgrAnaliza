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


int omegaTest(){
    TH1::SetDefaultSumw2();

    TFile* f = new TFile("mediumDataOmegaAR2.root","open");
    TDirectory *pbarOmega = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended"); //antyproton omega
    TDirectory *V0ChildPos = (TDirectory*)pbarOmega->Get("V0Child_pos");
    TDirectory *V0ChildNeg = (TDirectory*)pbarOmega->Get("V0Child_neg");
    TDirectory *omegaBachelor = (TDirectory*)pbarOmega->Get("hBachelor");
    TDirectory *omegaCascade = (TDirectory*)pbarOmega->Get("Cascade");

    TDirectory *pOmegaBar = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_antyomega"); //proton antyomega
    TDirectory *V0ChildPosAO = (TDirectory*)pOmegaBar->Get("V0Child_pos");
    TDirectory *V0ChildNegAO = (TDirectory*)pOmegaBar->Get("V0Child_neg");
    TDirectory *antyomegaBachelor = (TDirectory*)pOmegaBar->Get("hBachelor");
    TDirectory *antyomegaCascade = (TDirectory*)pOmegaBar->Get("Cascade");


    TH2F *hV0ChildPosPt = (TH2F*)V0ChildPos->Get("hPt");
    TH2F *hV0ChildNegPt = (TH2F*)V0ChildNeg->Get("hPt");
    TH2F *hOmegaBachelorPt = (TH2F*)omegaBachelor->Get("hPt");
    TH2F *hOmegaBachelorCharge = (TH2F*)omegaBachelor->Get("hCharge");
    TH2F *hOmegaCascadePt = (TH2F*)omegaCascade->Get("hPt");

    TH2F *hV0ChildPosPtAO = (TH2F*)V0ChildPosAO->Get("hPt");
    TH2F *hV0ChildNegPtAO = (TH2F*)V0ChildNegAO->Get("hPt");
    TH2F *hAntyOmegaBachelorPt = (TH2F*)antyomegaBachelor->Get("hPt");
    TH2F *hAntyOmegaBachelorCharge = (TH2F*)antyomegaBachelor->Get("hCharge");
    TH2F *hAntyOmegaCascadePt = (TH2F*)antyomegaCascade->Get("hPt");

    TCanvas *c1 = new TCanvas("pbarOmegaChildren","pbarOmegaChildren",800,600);
    c1->Divide(2,2);
    c1->cd(1);
    hV0ChildPosPt->Draw("COLZ");
    c1->cd(2);
    hV0ChildNegPt->Draw("COLZ");   
    c1->cd(3);
    hOmegaBachelorPt->Draw("COLZ");
    c1->cd(4);
    hOmegaBachelorCharge->Draw("COLZ");
    TCanvas *c2 = new TCanvas("pOmegaBarChildren","pOmegaBarChildren",800,600);
    c2->Divide(2,2);
    c2->cd(1);
    hV0ChildPosPtAO->Draw("COLZ");
    c2->cd(2);
    hV0ChildNegPtAO->Draw("COLZ");   
    c2->cd(3);
    hAntyOmegaBachelorPt->Draw("COLZ");
    c2->cd(4);
    hAntyOmegaBachelorCharge->Draw("COLZ");

    TCanvas *c3 = new TCanvas("cascade pt","cascade pt",800,600);
    c3->Divide(1,2);
    c3->cd(1);
    hAntyOmegaCascadePt->Draw("COLZ");
    c3->cd(2);
    hOmegaCascadePt->Draw("COLZ");
    

    return 0;

}