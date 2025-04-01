#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

TF1* performFit(TH1* hist, double* amp, double* mean, double* std, double start, double end){
    TF1* fitFunc = new TF1("fitFunc", "gaus", start, end);
    fitFunc->SetParameters(*amp,*mean,*std);

    hist->Fit(fitFunc,"R");
    *amp = fitFunc->GetParameter(0);
    *mean = fitFunc->GetParameter(1);
    *std = fitFunc->GetParameter(2);

    return fitFunc;
}



void fitPeaks(){
    // Enable batch mode (disables all pop-up windows)
    gROOT->SetBatch(kTRUE);

    // Set style to include fit parameters in the stats box
    gStyle->SetOptStat(1111);  // Show entries, mean, and RMS
    gStyle->SetOptFit(1111);   // Show fit results (Chi2, NDF, parameters)

    const char* rootFile = "./fitData/output_250221_Hole2_HV3600.root";
    const char* histName = "Cathode/charge/chargeL3";

    // Opening the ROOT file
    TFile* file = new TFile(rootFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return;
    }

    // Getting the histogram
    TH1* hist = (TH1*)file->Get(histName);
    if (!hist) {
        std::cerr << "Error: Could not retrieve histogram." << std::endl;
        file->Close();
        return;
    }

    hist->Rebin(128);
    
    // Get the inital fit parameters
    double max = hist->GetMaximum();
    double maxBin = hist->GetMaximumBin();
    double maxBinCenter = hist->GetBinCenter(maxBin);
    double std = 250.0; // Decent guess for the standard deviation

    double fitRange = 500;
    double min0 = maxBinCenter - fitRange;
    double max0 = maxBinCenter + fitRange;

    // Perform the fit
    TF1* fit = performFit(hist, &max, &maxBinCenter, &std, min0, max0);

    // Drawing the histogram
    TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
    hist->Draw();
    fit->Draw("same");

    // Save and close the file
    c->SaveAs("./fitData/complexFit/complexFitPeaks_rbf128.png");

    delete fit;
    delete c;

    file->Close();
    delete file;
}