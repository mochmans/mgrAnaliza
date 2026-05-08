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
#include "histSave.cc"

void cutsCompare(){


    TH1::SetDefaultSumw2();
    TString rootFileNameOld = "datasets/outputNew.root";
    TFile* fOld = new TFile(rootFileNameOld,"open");

    TString rootFileNameNew = "merged_refCuts.root";
    TFile* fNew = new TFile(rootFileNameNew,"open");

    double markerSize = 0.85;


    Int_t rebinFactor = 1;

    std::vector<TH1F*> omegaCorFuncsOld = correlationFunctions(fNew,"femto-universe-pair-task-track-cascade-extended_antyomega",rebinFactor);
    std::vector<TH1F*> antiOmegaCorFuncsOld = correlationFunctions(fOld,"femto-universe-pair-task-track-cascade-extended_antyomegaToF",rebinFactor);


    std::vector<TH1F*> omegaCorFuncsNew = correlationFunctions(fNew,"femto-universe-pair-task-track-cascade-extended",rebinFactor);
    std::vector<TH1F*> antiOmegaCorFuncsNew = correlationFunctions(fNew,"femto-universe-pair-task-track-cascade-extended_antyomega",rebinFactor);

    // std::vector<TH1F*> omegaCorFuncsNewToF = correlationFunctions(fNew,"femto-universe-pair-task-track-cascade-extended_allToF",rebinFactor);
    // std::vector<TH1F*> antiOmegaCorFuncsNewToF = correlationFunctions(fNew,"femto-universe-pair-task-track-cascade-extended_antyomegaToF",rebinFactor);


    TDirectory *DirNewProtons = (TDirectory*)fNew->Get("femto-universe-pair-task-track-cascade-extended_antyomega");
    TDirectory *DirNewAntiprotons = (TDirectory*)fNew->Get("femto-universe-pair-task-track-cascade-extended");


    TDirectory *sameEventProtons = (TDirectory*)DirNewProtons->Get("SameEvent");
    TDirectory *mixedEventProtons = (TDirectory*)DirNewProtons->Get("MixedEvent");


    TH2F *hAllCentrSameEventProtons = (TH2F*)sameEventProtons->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEventProtons = (TH2F*)mixedEventProtons->Get("relPairkstarMult");

    TH1D *hSameEventRatioProtons = hAllCentrSameEventProtons->ProjectionX("hSameEventRatioProtons",2,9);
    TH1D *hMixedEventRatioProtons = hAllCentrMixedEventProtons->ProjectionX("hMixedEventRatioProtons",2,9);

    TH1D *hkstarSameEventProtons0_10 = hAllCentrSameEventProtons->ProjectionX("hkstarSameEventProtons0_10",2,2);
    TH1D *hkstarMixedEventProtons0_10 = hAllCentrMixedEventProtons->ProjectionX("hkstarMixedEventProtons0_10",2,2);

    TH1D *hkstarSameEventProtons10_20 = hAllCentrSameEventProtons->ProjectionX("hkstarSameEventProtons10_20",3,3);
    TH1D *hkstarMixedEventProtons10_20 = hAllCentrMixedEventProtons->ProjectionX("hkstarMixedEventProtons10_20",3,3);

    TH1D *hkstarSameEventProtons20_30 = hAllCentrSameEventProtons->ProjectionX("hkstarSameEventProtons20_30",4,4);
    TH1D *hkstarMixedEventProtons20_30 = hAllCentrMixedEventProtons->ProjectionX("hkstarMixedEventProtons20_30",4,4);

    TH1D *hkstarSameEventProtons30_50 = hAllCentrSameEventProtons->ProjectionX("hkstarSameEventProtons30_50",5,6);
    TH1D *hkstarMixedEventProtons30_50 = hAllCentrMixedEventProtons->ProjectionX("hkstarMixedEventProtons30_50",5,6);

    TH1D *hkstarSameEventProtons50_100 = hAllCentrSameEventProtons->ProjectionX("hkstarSameEventProtons50_100",7,9);
    TH1D *hkstarMixedEventProtons50_100 = hAllCentrMixedEventProtons->ProjectionX("hkstarMixedEventProtons50_100",7,9);



    // ////////////////////////

    TDirectory *sameEventAntiprotons = (TDirectory*)DirNewAntiprotons->Get("SameEvent");
    TDirectory *mixedEventAntiprotons = (TDirectory*)DirNewAntiprotons->Get("MixedEvent");


    TH2F *hAllCentrSameEventAntiprotons = (TH2F*)sameEventAntiprotons->Get("relPairkstarMult");
    TH2F *hAllCentrMixedEventAntiprotons = (TH2F*)mixedEventAntiprotons->Get("relPairkstarMult");

    TH1D *hSameEventRatioAntiprotons = hAllCentrSameEventAntiprotons->ProjectionX("hSameEventRatioAntiprotons",2,9);
    TH1D *hMixedEventRatioAntiprotons = hAllCentrMixedEventAntiprotons->ProjectionX("hMixedEventRatioAntiprotons",2,9);

    TH1D *hkstarSameEventAntiprotons0_10 = hAllCentrSameEventAntiprotons->ProjectionX("hkstarSameEventAntiprotons0_10",2,2);
    TH1D *hkstarMixedEventAntiprotons0_10 = hAllCentrMixedEventAntiprotons->ProjectionX("hkstarMixedEventAntiprotons0_10",2,2);

    TH1D *hkstarSameEventAntiprotons10_20 = hAllCentrSameEventAntiprotons->ProjectionX("hkstarSameEventAntiprotons10_20",3,3);
    TH1D *hkstarMixedEventAntiprotons10_20 = hAllCentrMixedEventAntiprotons->ProjectionX("hkstarMixedEventAntiprotons10_20",3,3);

    TH1D *hkstarSameEventAntiprotons20_30 = hAllCentrSameEventAntiprotons->ProjectionX("hkstarSameEventAntiprotons20_30",4,4);
    TH1D *hkstarMixedEventAntiprotons20_30 = hAllCentrMixedEventAntiprotons->ProjectionX("hkstarMixedEventAntiprotons20_30",4,4);

    TH1D *hkstarSameEventAntiprotons30_50 = hAllCentrSameEventAntiprotons->ProjectionX("hkstarSameEventAntiprotons30_50",5,6);
    TH1D *hkstarMixedEventAntiprotons30_50 = hAllCentrMixedEventAntiprotons->ProjectionX("hkstarMixedEventAntiprotons30_50",5,6);

    TH1D *hkstarSameEventAntiprotons50_100 = hAllCentrSameEventAntiprotons->ProjectionX("hkstarSameEventAntiprotons50_100",7,9);
    TH1D *hkstarMixedEventAntiprotons50_100 = hAllCentrMixedEventAntiprotons->ProjectionX("hkstarMixedEventAntiprotons50_100",7,9);
    
    TH1F* hCControlProt = createCFHistogram((TH1F*)(hkstarSameEventProtons0_10), (TH1F*)(hkstarMixedEventProtons0_10), "hCControlProt", rebinFactor);   
    TH1F* hCControlAntiProt = createCFHistogram((TH1F*)(hkstarSameEventAntiprotons0_10), (TH1F*)(hkstarMixedEventAntiprotons0_10), "hCControlAntiProt", rebinFactor);
    std::cout << "Control histograms created" << std::endl;

    TH1F* hSameEventRatioMerged = (TH1F*)hSameEventRatioProtons->Clone("hSameEventRatioMerged");
    hSameEventRatioMerged->Add(hSameEventRatioAntiprotons);
    TH1F* hMixedEventRatioMerged = (TH1F*)hMixedEventRatioProtons->Clone("hMixedEventRatioMerged");
    hMixedEventRatioMerged->Add(hMixedEventRatioAntiprotons);

    TH1F* hkstarSameEventMerged0_10 = (TH1F*)hkstarSameEventProtons0_10->Clone("hkstarSameEventMerged0_10");
    hkstarSameEventMerged0_10->Add(hkstarSameEventAntiprotons0_10);
    TH1F* hkstarMixedEventMerged0_10 = (TH1F*)hkstarMixedEventProtons0_10->Clone("hkstarMixedEventMerged0_10");
    hkstarMixedEventMerged0_10->Add(hkstarMixedEventAntiprotons0_10);
    TH1F* hkstarSameEventMerged10_20 = (TH1F*)hkstarSameEventProtons10_20->Clone("hkstarSameEventMerged10_20");
    hkstarSameEventMerged10_20->Add(hkstarSameEventAntiprotons10_20);
    TH1F* hkstarMixedEventMerged10_20 = (TH1F*)hkstarMixedEventProtons10_20->Clone("hkstarMixedEventMerged10_20");
    hkstarMixedEventMerged10_20->Add(hkstarMixedEventAntiprotons10_20);
    TH1F* hkstarSameEventMerged20_30 = (TH1F*)hkstarSameEventProtons20_30->Clone("hkstarSameEventMerged20_30");
    hkstarSameEventMerged20_30->Add(hkstarSameEventAntiprotons20_30);
    TH1F* hkstarMixedEventMerged20_30 = (TH1F*)hkstarMixedEventProtons20_30->Clone("hkstarMixedEventMerged20_30");
    hkstarMixedEventMerged20_30->Add(hkstarMixedEventAntiprotons20_30);
    TH1F* hkstarSameEventMerged30_50 = (TH1F*)hkstarSameEventProtons30_50->Clone("hkstarSameEventMerged30_50");
    hkstarSameEventMerged30_50->Add(hkstarSameEventAntiprotons30_50);
    TH1F* hkstarMixedEventMerged30_50 = (TH1F*)hkstarMixedEventProtons30_50->Clone("hkstarMixedEventMerged30_50");
    hkstarMixedEventMerged30_50->Add(hkstarMixedEventAntiprotons30_50);
    TH1F* hkstarSameEventMerged50_100 = (TH1F*)hkstarSameEventProtons50_100->Clone("hkstarSameEventMerged50_100");
    hkstarSameEventMerged50_100->Add(hkstarSameEventAntiprotons50_100);
    TH1F* hkstarMixedEventMerged50_100 = (TH1F*)hkstarMixedEventProtons50_100->Clone("hkstarMixedEventMerged50_100");
    hkstarMixedEventMerged50_100->Add(hkstarMixedEventAntiprotons50_100);

    TH1F* hCFall = createCFHistogram((TH1F*)hSameEventRatioMerged, (TH1F*)hMixedEventRatioMerged, "hCFun", rebinFactor);   //antyproton omega
    TH1F* hCF0_10 = createCFHistogram((TH1F*)hkstarSameEventMerged0_10, (TH1F*)hkstarMixedEventMerged0_10, "hCF0_10", rebinFactor);
    TH1F* hCF10_20 = createCFHistogram((TH1F*)hkstarSameEventMerged10_20, (TH1F*)hkstarMixedEventMerged10_20, "hCF10_20", rebinFactor);
    TH1F* hCF20_30 = createCFHistogram((TH1F*)hkstarSameEventMerged20_30, (TH1F*)hkstarMixedEventMerged20_30, "hCF20_30", rebinFactor);
    TH1F* hCF30_50 = createCFHistogram((TH1F*)hkstarSameEventMerged30_50, (TH1F*)hkstarMixedEventMerged30_50, "hCF30_50", rebinFactor);
    TH1F* hCF50_100 = createCFHistogram((TH1F*)hkstarSameEventMerged50_100, (TH1F*)hkstarMixedEventMerged50_100, "hCF50_100", rebinFactor);

    TCanvas *cMerged = new TCanvas("cMerged","cMerged",800,600);
    cMerged->Divide(2,3);
    cMerged->cd(1);
    hCFall->Draw();
    hCFall->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCFall, kRed, kFullSquare, markerSize);
    cMerged->cd(2);
    hCF0_10->Draw();
    hCF0_10->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCF0_10, kRed, kFullSquare, markerSize);
    cMerged->cd(3);
    hCF10_20->Draw();
    hCF10_20->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCF10_20, kRed, kFullSquare, markerSize);
    cMerged->cd(4);
    hCF20_30->Draw();
    hCF20_30->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCF20_30, kRed, kFullSquare, markerSize);
    cMerged->cd(5);
    hCF30_50->Draw();
    hCF30_50->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCF30_50, kRed, kFullSquare, markerSize);
    cMerged->cd(6);
    hCF50_100->Draw();
    hCF50_100->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCF50_100, kRed, kFullSquare, markerSize);

    TCanvas *cControl = new TCanvas("cControl","cControl",800,600);
    cControl->Divide(1,2);
    cControl->cd(1);
    hCControlProt->Draw();
    hCControlProt->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCControlProt, kRed, kFullSquare, markerSize);
    cControl->cd(2);
    hCControlAntiProt->Draw();
    hCControlAntiProt->GetXaxis()->SetRangeUser(0,0.2);
    configureHistogram(hCControlAntiProt, kRed, kFullSquare, markerSize);

    //tof versions


    
    TDirectory *DirOldAntiprotons = (TDirectory*)fOld->Get("femto-universe-pair-task-track-cascade-extended_allToF");
    TDirectory *DirOldProtons = (TDirectory*)fOld->Get("femto-universe-pair-task-track-cascade-extended_antyomegaToF");
    // TDirectory *DirNewProtonsToF = (TDirectory*)fNew->Get("femto-universe-pair-task-track-cascade-extended_antyomegaToF");
    // TDirectory *DirNewAntiprotonsToF = (TDirectory*)fNew->Get("femto-universe-pair-task-track-cascade-extended_allToF");


    TString directoryNameProtons = "Tracks_pos";
    TString directoryNameAntiprotons = "Tracks_neg";

    TDirectory *oldAntiprotons = (TDirectory*)DirOldAntiprotons->Get(directoryNameAntiprotons);
    TDirectory *oldProtons = (TDirectory*)DirOldProtons->Get(directoryNameProtons);
    // TDirectory *newProtonsToF = (TDirectory*)DirNewProtonsToF->Get(directoryNameProtons);
    // TDirectory *newAntiprotonsToF = (TDirectory*)DirNewAntiprotonsToF->Get(directoryNameAntiprotons);
    TDirectory *newProtons = (TDirectory*)DirNewProtons->Get(directoryNameProtons);
    TDirectory *newAntiprotons = (TDirectory*)DirNewAntiprotons->Get(directoryNameAntiprotons);

    
    TH1D *hOldProtonsPt = (TH1D*)oldProtons->Get("hPt");
    TH1D *hOldAntiprotonsPt = (TH1D*)oldAntiprotons->Get("hPt");
    TH1D *hNewProtonsPt = (TH1D*)newProtons->Get("hPt");
    TH1D *hNewAntiprotonsPt = (TH1D*)newAntiprotons->Get("hPt");
    // TH1D *hNewProtonsPtToF = (TH1D*)newProtonsToF->Get("hPt");
    // TH1D *hNewAntiprotonsPtToF = (TH1D*)newAntiprotonsToF->Get("hPt");  

    TH2D *hOldProtonsnSigmaTPC = (TH2D*)oldProtons->Get("nSigmaTPC");
    TH2D *hOldAntiprotonsnSigmaTPC = (TH2D*)oldAntiprotons->Get("nSigmaTPC");
    TH2D *hNewProtonsnSigmaTPC = (TH2D*)newProtons->Get("nSigmaTPC");
    TH2D *hNewAntiprotonsnSigmaTPC = (TH2D*)newAntiprotons->Get("nSigmaTPC");
    // TH2D *hNewProtonsnSigmaTPCToF = (TH2D*)newProtonsToF->Get("nSigmaTPC");
    // TH2D *hNewAntiprotonsnSigmaTPCToF = (TH2D*)newAntiprotonsToF->Get("nSigmaTPC");  



    TDirectory *oldSameEventpAo = (TDirectory*)DirOldProtons->Get("SameEvent");
    TDirectory *oldMixedEventpAo = (TDirectory*)DirOldProtons->Get("MixedEvent");
    TDirectory *oldSameEventApo = (TDirectory*)DirOldAntiprotons->Get("SameEvent");
    TDirectory *oldMixedEventApo = (TDirectory*)DirOldAntiprotons->Get("MixedEvent");

    TH1D *hOldPairDistSameEventpAo = (TH1D*)oldSameEventpAo->Get("relPairDist");
    TH1D *hOldPairDistMixedEventpAo = (TH1D*)oldMixedEventpAo->Get("relPairDist");
    TH1D *hOldPairDistSameEventApo = (TH1D*)oldSameEventApo->Get("relPairDist");
    TH1D *hOldPairDistMixedEventApo = (TH1D*)oldMixedEventApo->Get("relPairDist");

    
    TDirectory *newSameEventpAo = (TDirectory*)DirNewProtons->Get("SameEvent");
    TDirectory *newMixedEventpAo = (TDirectory*)DirNewProtons->Get("MixedEvent");
    TDirectory *newSameEventApo = (TDirectory*)DirNewAntiprotons->Get("SameEvent");
    TDirectory *newMixedEventApo = (TDirectory*)DirNewAntiprotons->Get("MixedEvent");

    TH1D *hNewPairDistSameEventpAo = (TH1D*)newSameEventpAo->Get("relPairDist");
    TH1D *hNewPairDistMixedEventpAo = (TH1D*)newMixedEventpAo->Get("relPairDist");
    TH1D *hNewPairDistSameEventApo = (TH1D*)newSameEventApo->Get("relPairDist");
    TH1D *hNewPairDistMixedEventApo = (TH1D*)newMixedEventApo->Get("relPairDist");

    
    // TDirectory *newSameEventpAoToF = (TDirectory*)DirNewProtonsToF->Get("SameEvent");
    // TDirectory *newMixedEventpAoToF = (TDirectory*)DirNewProtonsToF->Get("MixedEvent");
    // TDirectory *newSameEventApoToF = (TDirectory*)DirNewAntiprotonsToF->Get("SameEvent");
    // TDirectory *newMixedEventApoToF = (TDirectory*)DirNewAntiprotonsToF->Get("MixedEvent");

    // TH1D *hNewPairDistSameEventpAoToF = (TH1D*)newSameEventpAoToF->Get("relPairDist");
    // TH1D *hNewPairDistMixedEventpAoToF = (TH1D*)newMixedEventpAoToF->Get("relPairDist");
    // TH1D *hNewPairDistSameEventApoToF = (TH1D*)newSameEventApoToF->Get("relPairDist");
    // TH1D *hNewPairDistMixedEventApoToF = (TH1D*)newMixedEventApoToF->Get("relPairDist");

    

    TCanvas *c12 = new TCanvas("Old cfs","Old cfs",800,600);
    c12->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(omegaCorFuncsOld[i], kBlue, kFullCircle, markerSize);
        configureHistogram(antiOmegaCorFuncsOld[i], kRed, kFullSquare, markerSize);
        omegaCorFuncsOld[i]->GetXaxis()->SetRangeUser(0,0.2);
        c12->cd(i+1);
        omegaCorFuncsOld[i]->Draw();
        antiOmegaCorFuncsOld[i]->Draw("SAME");
    }
    TLegend *legendOldCf = new TLegend();
    legendOldCf->AddEntry(omegaCorFuncsOld[0],"Omega ","p");
    legendOldCf->AddEntry(antiOmegaCorFuncsOld[0],"Antyomega ","p");
    c12->cd(1);
    legendOldCf->Draw();

    TCanvas *c13 = new TCanvas("New cfs","new cf",800,600);
    c13->Divide(2,3);
    for(Int_t i = 0; i<6; i++){
        configureHistogram(omegaCorFuncsNew[i], kRed, kFullSquare, markerSize);
        configureHistogram(omegaCorFuncsOld[i], kBlue, kFullCircle, markerSize);
        omegaCorFuncsNew[i]->GetXaxis()->SetRangeUser(0,0.2);
        c13->cd(i+1);
        omegaCorFuncsNew[i]->Draw();
        omegaCorFuncsOld[i]->Draw("SAME");  
    }
    TLegend *legendNewCf = new TLegend();
    legendNewCf->AddEntry(omegaCorFuncsNew[0],"New  ","p");
    legendNewCf->AddEntry(omegaCorFuncsOld[0],"Old ","p");
    c13->cd(1);
    legendNewCf->Draw();

    // TCanvas *c14 = new TCanvas("New ToF cfs","new ToF cf",800,600);
    // c14->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(omegaCorFuncsNewToF[i], kBlue, kFullCircle, markerSize);
    //     configureHistogram(antiOmegaCorFuncsNewToF[i], kRed, kFullSquare, markerSize);
    //     omegaCorFuncsNewToF[i]->GetXaxis()->SetRangeUser(0,0.2);
    //     c14->cd(i+1);
    //     omegaCorFuncsNewToF[i]->Draw();
    //     antiOmegaCorFuncsNewToF[i]->Draw("SAME");
    // }
    // TLegend *legendNewToFCf = new TLegend();
    // legendNewToFCf->AddEntry(omegaCorFuncsNewToF[0],"Omega antiproton ToF ","p");
    // legendNewToFCf->AddEntry(antiOmegaCorFuncsNewToF[0],"Antyomega proton ToF ","p");
    // c14->cd(1);
    // legendNewToFCf->Draw(); 

    // TCanvas *c15 = new TCanvas("ToF 0.5 vs ToF 0.75");
    // c15->Divide(2,3);
    // for(Int_t i = 0; i<6; i++){
    //     configureHistogram(omegaCorFuncsNew[i], kRed, kFullSquare, markerSize);
    //     configureHistogram(omegaCorFuncsNewToF[i], kBlue, kFullCircle, markerSize);
    //     omegaCorFuncsNewToF[i]->GetXaxis()->SetRangeUser(0,0.2);
    //     c15->cd(i+1);
    //     omegaCorFuncsNewToF[i]->Draw();
    //     antiOmegaCorFuncsNewToF[i]->Draw("SAME");
    // }
    // TLegend *legendToFvsNoToF = new TLegend();
    // legendToFvsNoToF->AddEntry(omegaCorFuncsNewToF[0],"Omega antiproton ToF 0.5 ","p");
    // legendToFvsNoToF->AddEntry(omegaCorFuncsNew[0],"Omega proton ToF 0.75 ","p");
    // c15->cd(1);
    // legendToFvsNoToF->Draw();

    TCanvas *c21 = new TCanvas("K*","Old pair dist",800,600);
    c21->Divide(2,1);
    configureHistogram(hOldPairDistSameEventpAo, kBlue, kFullCircle, markerSize);
    configureHistogram(hOldPairDistMixedEventpAo, kBlue, kFullCircle, markerSize);

    configureHistogram(hOldPairDistSameEventApo, kRed, kFullSquare, markerSize);
    configureHistogram(hOldPairDistMixedEventApo, kRed, kFullSquare, markerSize);
    
    configureHistogram(hNewPairDistSameEventApo, kGreen+2, kFullCross, markerSize);
    configureHistogram(hNewPairDistMixedEventApo, kGreen+2, kFullCross, markerSize);

    configureHistogram(hNewPairDistSameEventpAo, kMagenta, kFullDiamond, markerSize);
    configureHistogram(hNewPairDistMixedEventpAo, kMagenta, kFullDiamond, markerSize);
    c21->cd(1);
    hNewPairDistSameEventpAo->Draw();
    hOldPairDistSameEventpAo->Draw("SAME");
    hOldPairDistSameEventApo->Draw("SAME");
    hNewPairDistSameEventApo->Draw("SAME");
    
    TLegend *legendOldPairDist = new TLegend();
    legendOldPairDist->AddEntry(hOldPairDistSameEventpAo," p-antiomega","p");
    legendOldPairDist->AddEntry(hNewPairDistSameEventpAo," p-antiomega new","p");
    legendOldPairDist->AddEntry(hOldPairDistSameEventApo," anti-p-Omega","p");
    legendOldPairDist->AddEntry(hNewPairDistSameEventApo," anti-p-Omega new","p");
    legendOldPairDist->Draw();
    c21->cd(2);
    hNewPairDistMixedEventpAo->Draw();
    hOldPairDistMixedEventpAo->Draw("SAME");
    hOldPairDistMixedEventApo->Draw("SAME");
    hNewPairDistMixedEventApo->Draw("SAME");


    // TCanvas *c22 = new TCanvas("K* ToF","New pair dist ToF",800,600);
    // c22->Divide(2,1);
    // configureHistogram(hNewPairDistSameEventpAoToF, kBlue, kFullDiamond, markerSize);
    // configureHistogram(hNewPairDistMixedEventpAoToF, kBlue, kFullDiamond, markerSize);
        

    // c22->cd(1);
    
    // hNewPairDistSameEventpAo->Draw();
    // hNewPairDistSameEventpAoToF->Draw("SAME");
    // TLegend *legendToFvsNoToFPairDistSame = new TLegend();
    // legendToFvsNoToFPairDistSame->AddEntry(hNewPairDistSameEventpAoToF," p - antiomega ToF 0.5 ","p");
    // legendToFvsNoToFPairDistSame->AddEntry(hNewPairDistSameEventpAo," p - antiomega ToF 0.75 ","p"); 
    // legendToFvsNoToFPairDistSame->Draw();   
    // c22->cd(2);
    
    // hNewPairDistMixedEventpAo->Draw();
    // hNewPairDistMixedEventpAoToF->Draw("SAME");

    // TCanvas *c23 = new TCanvas("nSigmaTPC","nSigmaTPC",800,600);
    // c23->Divide(2,3);
    // c23->cd(1);
    // hOldProtonsnSigmaTPC->Draw("COLZ");
    // hOldProtonsnSigmaTPC->GetXaxis()->SetRangeUser(0.5,4);
    // hOldProtonsnSigmaTPC->GetYaxis()->SetRangeUser(-3,3);
    // hOldProtonsnSigmaTPC->SetStats(0);
    // c23->cd(2);
    // hOldAntiprotonsnSigmaTPC->Draw("COLZ"); 
    // hOldAntiprotonsnSigmaTPC->GetXaxis()->SetRangeUser(0.5,4);
    // hOldAntiprotonsnSigmaTPC->GetYaxis()->SetRangeUser(-3,3);
    // hOldAntiprotonsnSigmaTPC->SetStats(0);

    // c23->cd(3);
    // hNewProtonsnSigmaTPC->Draw("COLZ");
    // hNewProtonsnSigmaTPC->GetXaxis()->SetRangeUser(0.5,4);
    // hNewProtonsnSigmaTPC->GetYaxis()->SetRangeUser(-3,3);
    //     hNewProtonsnSigmaTPC->SetStats(0);
    // c23->cd(4);
    // hNewAntiprotonsnSigmaTPC->Draw("COLZ");
    // hNewAntiprotonsnSigmaTPC->GetXaxis()->SetRangeUser(0.5,4);
    // hNewAntiprotonsnSigmaTPC->GetYaxis()->SetRangeUser(-3,3);
    // hNewAntiprotonsnSigmaTPC->SetStats(0);
    // c23->cd(5);
    // hNewProtonsnSigmaTPCToF->Draw("COLZ");
    // hNewProtonsnSigmaTPCToF->GetXaxis()->SetRangeUser(0.5,4);
    // hNewProtonsnSigmaTPCToF->GetYaxis()->SetRangeUser(-3,3);
    // hNewProtonsnSigmaTPCToF->SetStats(0);
    // c23->cd(6); 
    // hNewAntiprotonsnSigmaTPCToF->Draw("COLZ");
    // hNewAntiprotonsnSigmaTPCToF->GetXaxis()->SetRangeUser(0.5,4);
    // hNewAntiprotonsnSigmaTPCToF->GetYaxis()->SetRangeUser(-3,3);
    // hNewAntiprotonsnSigmaTPCToF->SetStats(0);

    TCanvas *c24 = new TCanvas("Pt","Pt",800,600);
    c24->Divide(2,1);

    configureHistogram(hOldProtonsPt, kBlue, kFullCircle, markerSize);
    configureHistogram(hOldAntiprotonsPt, kBlue, kFullSquare, markerSize);
    configureHistogram(hNewProtonsPt, kRed, kFullCircle, markerSize);
    configureHistogram(hNewAntiprotonsPt, kRed, kFullSquare, markerSize);
    // configureHistogram(hNewProtonsPtToF, kGreen+2, kFullDiamond, markerSize);
    // configureHistogram(hNewAntiprotonsPtToF, kMagenta, kFullCross, markerSize);

    c24->cd(1);
    hNewProtonsPt->Draw();
    hOldProtonsPt->Draw("SAME");
    TLegend *legendPtOld = new TLegend();
    legendPtOld->AddEntry(hOldProtonsPt,"(anti)Protons old","p");
    legendPtOld->AddEntry(hNewProtonsPt,"(anti)Protons new","p");
    legendPtOld->Draw();
    c24->cd(2);
    hNewAntiprotonsPt->Draw();
    hOldAntiprotonsPt->Draw("SAME");
    // c24->cd(3);
    // hNewProtonsPtToF->Draw();
    // hNewAntiprotonsPtToF->Draw("SAME");

    // c24->cd(4);
    // hNewAntiprotonsPt->Draw();
    // // hNewAntiprotonsPtToF->Draw("SAME");
    // hNewProtonsPt->Draw("SAME");
    // // hNewProtonsPtToF->Draw("SAME");
    // TLegend *legendPtNew = new TLegend();
    // // legendPtNew->AddEntry(hNewProtonsPt,"Protons ToF 0.75","p");
    // // legendPtNew->AddEntry(hNewProtonsPtToF,"Protons ToF 0.5","p");
    // // legendPtNew->AddEntry(hNewAntiprotonsPt,"Antiprotons ToF 0.75","p");
    // // legendPtNew->AddEntry(hNewAntiprotonsPtToF,"Antiprotons ToF 0.5","p");
    // legendPtNew->Draw();

    

}
    
    