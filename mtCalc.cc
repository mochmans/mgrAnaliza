#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TDirectory.h"

TH2F* getPtKstarDist(TString filename, int part = 1, bool antiPair = 0, bool sameEvent = 1)
{
    TFile* f = new TFile(filename, "open");
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir;
    TDirectory *event;

    if(!antiPair)
        Dir = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended");
    else 
        Dir = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_antyomega");


    if(sameEvent)
        event = (TDirectory*)Dir->Get("SameEvent");
    else 
        event = (TDirectory*)Dir->Get("MixedEvent");

    if(part != 1 && part != 2)
    {
        std::cerr << "Invalid part number. Use 1 for Proton and 2 for Omega. Using proton." << std::endl;
        part = 1;
    }

    TH2F *hPtKstar = (TH2F*)event->Get(Form("kstarPtPart%d", part));

    return hPtKstar;

}

TH2F* getKtKstarDist(TString filename, bool antiPair = 0, bool sameEvent = 1){
    TFile* f = new TFile(filename, "open");
    
    TH1::SetDefaultSumw2();

    TDirectory *Dir;
    TDirectory *event;

    if(!antiPair)
        Dir = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended");
    else 
        Dir = (TDirectory*)f->Get("femto-universe-pair-task-track-cascade-extended_antyomega");
    if(sameEvent)
        event = (TDirectory*)Dir->Get("SameEvent");
    else 
        event = (TDirectory*)Dir->Get("MixedEvent");
    
    TH2F *hKtKstar = (TH2F*)event->Get("relPairkstarkT");

    return hKtKstar;
}

TH2F* addNewDists(TH2F* dist1, TH2F* dist2)
{
    if(!dist1 || !dist2)
    {
        std::cerr << "One of the distributions is null. Cannot add." << std::endl;
        return nullptr;
    }

    TH2F* result = (TH2F*)dist1->Clone("added dist");
    result->Add(dist2);
    return result;
}

void addDists(TH2F* dist1, TH2F* dist2)
{
    if(!dist1 || !dist2)
    {
        std::cerr << "One of the distributions is null. Cannot add." << std::endl;
        return;
    }

    dist1->Add(dist2);
}

void calculateAvgKt(){
    TH2F* hkTkStarzzh = getKtKstarDist("zzh_medium.root", false, true);
    TH2F* hkTkStarzzm = getKtKstarDist("zzm_medium.root", false, true);
    TH2F* hkTkStarar = getKtKstarDist("mediumDataOmegaARfull.root", false, true);
    TH2F* hkTkStarzzo = getKtKstarDist("zzo_medium.root", false, true);

    TH2F* hkTkStarAntiPairzzh = getKtKstarDist("zzh_medium.root", true, true);
    TH2F* hkTkStarAntiPairzzm = getKtKstarDist("zzm_medium.root", true, true);
    TH2F* hkTkStarAntiPairar = getKtKstarDist("mediumDataOmegaARfull.root", true, true);
    TH2F* hkTkStarAntiPairzzo = getKtKstarDist("zzo_medium.root", true, true);

    TH2F* hkTkStarPairSum = addNewDists(addNewDists(hkTkStarzzh, hkTkStarzzm), addNewDists(hkTkStarar, hkTkStarzzo));
    TH2F* hkTkStarAntiPairSum = addNewDists(addNewDists(hkTkStarAntiPairzzh, hkTkStarAntiPairzzm), addNewDists(hkTkStarAntiPairar, hkTkStarAntiPairzzo));

    TH2F* hkTkStarPairSumMerged = addNewDists(hkTkStarPairSum, hkTkStarAntiPairSum);

    double kStarMax = 0.15;

    TH1D* hAvgKtPair = hkTkStarPairSum->ProjectionY("hAvgKtPair", 1, hkTkStarPairSum->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgKtAntiPair = hkTkStarAntiPairSum->ProjectionY("hAvgKtAntiPair", 1, hkTkStarAntiPairSum->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgKtPairMerged = hkTkStarPairSumMerged->ProjectionY("hAvgKtPairMerged", 1, hkTkStarPairSumMerged->GetXaxis()->FindBin(kStarMax));

    TCanvas* c = new TCanvas("kT in k* range","kT in k* range", 800, 600);
    c->Divide(1,2);
    c->cd(1);
    hAvgKtPair->SetMarkerStyle(kFullSquare);
    hAvgKtPair->SetMarkerColor(kBlue);
    hAvgKtPair->SetMarkerSize(0.5);
    hAvgKtPair->SetStats(1);
    hAvgKtPair->Draw();    
    c->cd(2);
    hAvgKtAntiPair->SetMarkerStyle(kFullCircle);
    hAvgKtAntiPair->SetMarkerColor(kRed);
    hAvgKtAntiPair->SetMarkerSize(0.5);
    hAvgKtAntiPair->SetStats(1);
    hAvgKtAntiPair->Draw();
    TCanvas* c2 = new TCanvas("kT in k* range merged","kT in k* range merged", 800, 600);
    c2->cd(1);
    hAvgKtPairMerged->SetMarkerStyle(kFullSquare);
    hAvgKtPairMerged->SetMarkerColor(kBlue);
    hAvgKtPairMerged->SetMarkerSize(0.5);
    hAvgKtPairMerged->SetStats(1);
    hAvgKtPairMerged->Draw();

    
}


void calculateAvgPt()
{
    TH2F* hPtKstarAntiProtonzzh = getPtKstarDist("zzh_medium.root", 1, false, true);
    TH2F* hPtKstarAntiProtonzzm = getPtKstarDist("zzm_medium.root", 1, false, true);
    TH2F* hPtKstarAntiProtonar = getPtKstarDist("mediumDataOmegaARfull.root", 1, false, true);
    TH2F* hPtKstarAntiProtonzzo = getPtKstarDist("zzo_medium.root", 1, false, true);

    TH2F* hPtKstarOmegazzh = getPtKstarDist("zzh_medium.root", 2, false, true);
    TH2F* hPtKstarOmegazzm = getPtKstarDist("zzm_medium.root", 2, false, true);
    TH2F* hPtKstarOmegaar = getPtKstarDist("mediumDataOmegaARfull.root", 2, false, true);
    TH2F* hPtKstarOmegazzo = getPtKstarDist("zzo_medium.root", 2, false, true);

    TH2F* hPtKstarProtonzzh = getPtKstarDist("zzh_medium.root", 1, true, true);
    TH2F* hPtKstarProtonzzm = getPtKstarDist("zzm_medium.root", 1, true, true);
    TH2F* hPtKstarProtonar = getPtKstarDist("mediumDataOmegaARfull.root", 1, true, true);       
    TH2F* hPtKstarProtonzzo = getPtKstarDist("zzo_medium.root", 1, true, true);

    TH2F* hPtKstarAntiOmegazzh = getPtKstarDist("zzh_medium.root", 2, true, true);
    TH2F* hPtKstarAntiOmegazzm = getPtKstarDist("zzm_medium.root", 2, true, true);  
    TH2F* hPtKstarAntiOmegaar = getPtKstarDist("mediumDataOmegaARfull.root", 2, true, true);
    TH2F* hPtKstarAntiOmegazzo = getPtKstarDist("zzo_medium.root", 2, true, true);

    TH2F* hPtKstarAntiProtonSum = addNewDists(addNewDists(hPtKstarAntiProtonzzh, hPtKstarAntiProtonzzm), addNewDists(hPtKstarAntiProtonar, hPtKstarAntiProtonzzo));
    TH2F* hPtKstarOmegaSum = addNewDists(addNewDists(hPtKstarOmegazzh, hPtKstarOmegazzm), addNewDists(hPtKstarOmegaar, hPtKstarOmegazzo));

    TH2F* hPtKstarProtonSum = addNewDists(addNewDists(hPtKstarProtonzzh, hPtKstarProtonzzm), addNewDists(hPtKstarProtonar, hPtKstarProtonzzo));
    TH2F* hPtKstarAntiOmegaSum = addNewDists(addNewDists(hPtKstarAntiOmegazzh, hPtKstarAntiOmegazzm), addNewDists(hPtKstarAntiOmegaar, hPtKstarAntiOmegazzo));

    TH2F* hPtKstarProtonSumMerged = addNewDists(hPtKstarProtonSum, hPtKstarAntiProtonSum);
    TH2F* hPtKstarOmegaSumMerged = addNewDists(hPtKstarOmegaSum, hPtKstarAntiOmegaSum);

    double kStarMax = 0.15;

    TH1D* hAvgPtProton = hPtKstarProtonSum->ProjectionY("hAvgPtProton", 1, hPtKstarProtonSum->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgPtOmega = hPtKstarOmegaSum->ProjectionY("hAvgPtOmega", 1, hPtKstarOmegaSum->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgPtAntiProton = hPtKstarAntiProtonSum->ProjectionY("hAvgPtAntiProton", 1, hPtKstarAntiProtonSum->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgPtAntiOmega = hPtKstarAntiOmegaSum->ProjectionY("hAvgPtAntiOmega", 1, hPtKstarAntiOmegaSum->GetXaxis()->FindBin(kStarMax));

    TH1D* hAvgPtProtonSumMerged = hPtKstarProtonSumMerged->ProjectionY("hAvgPtProtonSumMerged", 1, hPtKstarProtonSumMerged->GetXaxis()->FindBin(kStarMax));
    TH1D* hAvgPtOmegaSumMerged = hPtKstarOmegaSumMerged->ProjectionY("hAvgPtOmegaSumMerged", 1, hPtKstarOmegaSumMerged->GetXaxis()->FindBin(kStarMax));

    TCanvas* c = new TCanvas("Average Pt", "Average Pt", 800, 600);
    c->Divide(2,2);
    c->cd(1);
    hAvgPtProton->SetMarkerColor(kBlue);
    hAvgPtProton->SetMarkerStyle(kFullCircle);
    hAvgPtProton->SetMarkerSize(0.5);
    hAvgPtProton->SetStats(1);
    hAvgPtProton->Draw();

    c->cd(2);
    hAvgPtOmega->SetMarkerColor(kRed);
    hAvgPtOmega->SetMarkerStyle(kFullSquare);
    hAvgPtOmega->SetMarkerSize(0.5);
    hAvgPtOmega->SetStats(1);
    hAvgPtOmega->Draw();

    c->cd(3);
    hAvgPtAntiProton->SetMarkerColor(kBlue);
    hAvgPtAntiProton->SetMarkerStyle(kFullCircle);
    hAvgPtAntiProton->SetMarkerSize(0.5);
    hAvgPtAntiProton->SetStats(1);
    hAvgPtAntiProton->Draw();

    c->cd(4);
    hAvgPtAntiOmega->SetMarkerColor(kRed);
    hAvgPtAntiOmega->SetMarkerStyle(kFullSquare);
    hAvgPtAntiOmega->SetMarkerSize(0.5);
    hAvgPtAntiOmega->SetStats(1);
    hAvgPtAntiOmega->Draw();

    TCanvas* c2 = new TCanvas("Merged Average Pt", "Merged Average Pt", 800, 600);
    c2->Divide(1,2);
    c2->cd(1);
    hAvgPtProtonSumMerged->SetMarkerColor(kBlue);
    hAvgPtProtonSumMerged->SetMarkerStyle(kFullCircle);
    hAvgPtProtonSumMerged->SetMarkerSize(0.5);
    hAvgPtProtonSumMerged->SetStats(1);
    hAvgPtProtonSumMerged->Draw();  

    c2->cd(2);
    hAvgPtOmegaSumMerged->SetMarkerColor(kRed);
    hAvgPtOmegaSumMerged->SetMarkerStyle(kFullSquare);
    hAvgPtOmegaSumMerged->SetMarkerSize(0.5);
    hAvgPtOmegaSumMerged->SetStats(1);
    hAvgPtOmegaSumMerged->Draw();   

}