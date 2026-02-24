#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include <stdio.h>
#include <iostream>
#include "../mergedCFs.cc"

void FitCorrMacro(TString filename, TString coulFileName,  Int_t centr = 0,Int_t n=1,const int Mbins=4){

    TCanvas *c1 = new TCanvas("c1", "Plots", 10, 10, 600, 500);
    

    std::vector<TH1F*> correctedCFs;
    HistogramGroups mergedHists = loadkStarHistograms(filename);

    std::cout<< mergedHists.omegaCor.size() << std::endl; // assuming all have the same number of histograms
    HistogramGroups mergedkStarMergedPairAntipair = kStarMergedPairAntipair(mergedHists);   //zmienia histogramy mergedHits

    // experimental function
    // TH1D* kstarSame = ((TH1D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended/SameEvent/relPairDist")));
    // TH1D* kstarMixed = ((TH1D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended/MixedEvent/relPairDist")));
    // TH1D* kstarSame2 = ((TH1D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/SameEvent/relPairDist")));
    // TH1D* kstarMixed2 = ((TH1D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/MixedEvent/relPairDist")));
    
    // kstarSame->Add(kstarSame2);
    // kstarMixed->Add(kstarMixed2);

    // zzh
    // TH2D* kstarSameMult = ((TH2D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended/SameEvent/relPairkstarMult")));
    // TH2D* kstarMixedMult = ((TH2D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended/MixedEvent/relPairkstarMult")));
    // TH2D* kstarSameMult2 = ((TH2D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/SameEvent/relPairkstarMult")));
    // TH2D* kstarMixedMult2 = ((TH2D*)ifile1->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/MixedEvent/relPairkstarMult")));
    // // zzl
    // TH2D* kstarSameMulti = ((TH2D*)ifile3->Get(Form("femto-universe-pair-task-track-cascade-extended/SameEvent/relPairkstarMult")));
    // TH2D* kstarMixedMulti = ((TH2D*)ifile3->Get(Form("femto-universe-pair-task-track-cascade-extended/MixedEvent/relPairkstarMult")));
    // TH2D* kstarSameMult2i = ((TH2D*)ifile3->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/SameEvent/relPairkstarMult")));
    // TH2D* kstarMixedMult2i = ((TH2D*)ifile3->Get(Form("femto-universe-pair-task-track-cascade-extended_p_axi/MixedEvent/relPairkstarMult")));
    
    // kstarSameMult->Add(kstarSameMult2);
    // kstarMixedMult->Add(kstarMixedMult2);

    // kstarSameMult->Add(kstarSameMulti);
    // kstarMixedMult->Add(kstarMixedMulti);
    // kstarSameMult->Add(kstarSameMult2i);
    // kstarMixedMult->Add(kstarMixedMult2i);

    // TH1D* kstarSame = ((TH1D*)ifile1->Get(Form("relPairDistSame")));
    // TH1D* kstarMixed = ((TH1D*)ifile1->Get(Form("relPairDistMixed")));

    // TH2D* kstarSameMult = ((TH2D*)ifile1->Get(Form("relPairkstarMultSame")));
    // TH2D* kstarMixedMult = ((TH2D*)ifile1->Get(Form("relPairkstarMultMixed")));

    TH1F* kstarSame = nullptr;
    TH1F* kstarMixed = nullptr;

    // TH1D* kstarSame = kstarSameMult->ProjectionX("same", 2, 2);
    // TH1D* kstarMixed = kstarMixedMult->ProjectionX("mixed", 2, 2);

    if(centr == 0){ //all centralities
        std::cout << "Using all centralities" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[0]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[1]->Clone("mixed"));
    }
    else if(centr == 1){  //0-10
        std::cout << "Using 0-10\% centrality" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[2]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[3]->Clone("mixed"));
    }
    else if(centr == 2){    //10-20
        std::cout << "Using 10-20\% centrality" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[4]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[5]->Clone("mixed"));
    }
    else if(centr == 3){    //20-30
        std::cout << "Using 20-30\% centrality" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[6]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[7]->Clone("mixed"));
    }
    else if(centr == 4){    //30-50
        std::cout << "Using 30-50\% centrality" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[8]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[9]->Clone("mixed"));
    }
    else if(centr == 5){    //50-100
        std::cout << "Using 50-100\% centrality" << std::endl;
        kstarSame = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[10]->Clone("same"));
        kstarMixed = ((TH1F*)mergedkStarMergedPairAntipair.omegaCor[11]->Clone("mixed"));
    }
    else{
        std::cerr << "Invalid centrality bin: " << centr << std::endl;
        return;
    }
    

    // TH1D* kstarSame2 = kstarSameMult->ProjectionX("same", 7, 7);
    // TH1D* kstarMixed2 = kstarMixedMult->ProjectionX("mixed", 7, 7);

    // kstarSame->Add(kstarSame2);
    // kstarMixed->Add(kstarMixed2);

    // kstarSame->Rebin(5);
    // kstarMixed->Rebin(5);

    TH1F* divideKstar = createCFHistogram(kstarSame, kstarMixed, "kstarCF", 5);

    // double scaleSame = kstarSame->Integral(kstarSame->FindBin(0.3), kstarSame->FindBin(0.5));
    // scaleSame = scaleSame/(kstarSame->FindBin(0.5)-kstarSame->FindBin(0.3));
    // double scaleMixed = kstarMixed->Integral(kstarMixed->FindBin(0.3), kstarMixed->FindBin(0.5));
    // scaleMixed = scaleMixed/(kstarMixed->FindBin(0.5)-kstarMixed->FindBin(0.3));

    // divideKstar->Scale(1/scaleSame);
    // kstarMixed->Scale(1/scaleMixed);

    // divideKstar->Divide(kstarMixed);

    // string rvalues[9] = {"6.00", "6.25", "6.50", "6.75", "7.00", "7.25", "7.50", "7.75", "8.00"};
    // double chimin = 999999;
    // string Rmin = "";

    
    TFile *ifile2 = new TFile(coulFileName);
    TGraph* coulFit = ((TGraph*)ifile2->Get("KCoulomb"));
    ifile2->Close();

    // cropping to reasonable fit
    double xmin = 0.0;
    double xmax = 0.13;

    std::vector<double> x_vals;
    std::vector<double> y_vals;

    double x1, y1;

    for (int i = 0; i < coulFit->GetN(); ++i) {
        coulFit->GetPoint(i, x1, y1);
        if (x1 >= xmin && x1 <= xmax) {
            x_vals.push_back(x1);
            y_vals.push_back(y1);
        }
    }
    // adding dummy points to converge at 1
    x_vals.push_back(0.3);
    y_vals.push_back(1.0);
    x_vals.push_back(0.4);
    y_vals.push_back(1.0);
    x_vals.push_back(0.5);
    y_vals.push_back(1.0);
    x_vals.push_back(0.6);
    y_vals.push_back(1.0);

    TGraph* coulFit_cropp = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());

    // spline extrapolation
    TSpline3 *spline = new TSpline3("spline", coulFit_cropp);
    TF1 *interpFunc = new TF1("interpFunc", [spline](double *x, double *) {
        return spline->Eval(x[0]);
    }, 0., 0.6, 0);

    // auto N = divideKstar->GetNbinsX();

    // change graph to histogram - OLD
    // auto dN = ((float)coulFit->GetN())/((float)divideKstar->GetNbinsX()); // number of points in TGraph / N of bins in data hist
    // TH1D* hcoulFit = new TH1D("hcoulFit", "HIST", 4*N, 0, 0.5);

    // for(float i=0; i < N; i+=dN) {
    //     //cout << "fit N:" << coulFit->GetN() << "; N:" << N << "; dN:" << dN << "; i:" << i << endl;
    //     double x,y;
    //     coulFit->GetPoint(i, x, y);
    //     //cout << y << endl;
    //     // if (y > 1.001){
    //     //     y=1;
    //     // }
    //     cout << y << endl;
    //     hcoulFit->Fill(x,y);
    // }

    // change to histogram
    // TH1D* hcoulFit = new TH1D("hcoulFit", "Histogram from TF1", N, 0.0, 0.6);
    // for (int i = 1; i <= N; i++) {
    //     double x = hcoulFit->GetBinCenter(i);
    //     double y = interpFunc->Eval(x);
    //     hcoulFit->SetBinContent(i, y);
    // }

    // poprawka na purity
    double p_purity = 0.98;
    double xi_purity = 0.8;
    double purity = 1;// p_purity*xi_purity; // pair purity

    TF1 *fscaled = new TF1("fscaled", [interpFunc](double *x, double *p) {
        return (interpFunc->Eval(x[0]) - 1.0) * p[0] + 1.0; 
    }, interpFunc->GetXmin(), interpFunc->GetXmax(), 1);

    fscaled->SetParameter(0, purity);


    divideKstar->SetTitle("Coulomb fit;k* (GeV/c);C(k*) ");
    divideKstar->GetXaxis()->SetLabelSize(0.045);
    divideKstar->GetXaxis()->SetTitleSize(0.05);
    divideKstar->GetYaxis()->SetLabelSize(0.045);
    divideKstar->GetYaxis()->SetTitleSize(0.05);
    divideKstar->GetYaxis()->SetRangeUser(0.3, 1.3);
    divideKstar->GetXaxis()->SetRangeUser(0., 0.5);
    divideKstar->SetMarkerColor(kGreen-3);
    divideKstar->SetMarkerStyle(8);
    divideKstar->SetMarkerSize(1);
    divideKstar->Draw("p");
    // hcoulFit->SetMarkerColor(kViolet);
    // hcoulFit->SetMarkerStyle(8);
    // hcoulFit->SetMarkerSize(0.6);
    // hcoulFit->SetLineColor(kViolet);
    // hcoulFit->SetLineWidth(2);
    // hcoulFit->Draw("p same");
    //interpFunc->SetLineColor(kBlue);
    //interpFunc->Draw("same");
    fscaled->Draw("same");

    auto Legend = new TLegend(0.6, 0.2, .8, .35);
    //Legend->AddEntry(interpFunc,"unscaled","L");
    Legend->AddEntry(divideKstar,"data","p");
    Legend->AddEntry(fscaled,"Coulomb fit","L");
    Legend->SetBorderSize(0);
    Legend->Draw();

    TLegend *text1 = new TLegend(0.65, 0.75, 0.8, 0.85);
    text1->SetTextSize(.06);
    if(centr == 0)
        text1->SetHeader("0-100\% centrality", "c");
    else if(centr == 1)
        text1->SetHeader("0-10\% centrality", "c");
    else if(centr == 2)
        text1->SetHeader("10-20\% centrality", "c");
    else if(centr == 3)
        text1->SetHeader("20-30\% centrality", "c");
    else if(centr == 4)
        text1->SetHeader("30-50\% centrality", "c");
    else if(centr == 5)
        text1->SetHeader("50-100\% centrality", "c");
    text1->SetBorderSize(0);
    text1->SetFillColor(0);
    text1->Draw("same");

    // TLegend *text22 = new TLegend(0.6, 0.35, 0.75, 0.5);
    // text22->SetTextSize(.06);
    // text22->SetTextFont(62);
    // text22->SetHeader("R = 0.01", "c");
    // text22->SetBorderSize(0);
    // text22->SetFillColor(0);
    // text22->Draw("same");

    // chi2
    double chi2 = 0;
    // chi2 = divideKstar->Chi2Test(hcoulFit, " CHI2");
    chi2 = divideKstar->Chisquare(fscaled,"");
    cout <<"chi2 = " << chi2 << endl;
    
}
