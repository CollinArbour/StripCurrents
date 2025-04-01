#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void fitPeaksComplex(){
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

    int rebinFactor = 256;
    hist->Rebin(rebinFactor);
    
    // Get the inital fit parameters
    double max = hist->GetMaximum();
    double max2 = max/2;
    double maxBin = hist->GetMaximumBin();
    double maxBinCenter = hist->GetBinCenter(maxBin);
    double std = 100.0; // Decent guess for the standard deviation
    double std2 = 300.0; // Decent guess for the standard deviation

    double fitRange = 800;
    double min0 = maxBinCenter - 1*fitRange;
    double max0 = maxBinCenter + 2*fitRange;

    /*
    // Fitting two gaussians
    // TF1* fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3)",1000,3000);
    TF1* fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3)",min0,max0);
    fitFunc->SetParameters(max,maxBinCenter,std,max2,maxBinCenter,std2);
    */

    /*
    // Fitting Crystal Ball Function
    TF1* fitFunc = new TF1("fitFunc","crystalball",500,max0);
    fitFunc->SetParameters(max,maxBinCenter,std,1500,2);
    */

    // Fitting Gaus + Landau
    TF1* fitFunc = new TF1("fitFunc","gaus(0) + landau(3)",min0,max0);
    fitFunc->SetParameters(max2,maxBinCenter,std2,max,maxBinCenter,std);

    /*
    // Fitting  Landau
    // VERY GOOD FIT
    TF1* fitFunc = new TF1("fitFunc","landau",min0,max0);
    fitFunc->SetParameters(max,maxBinCenter,std2);
    */
    
    /*
    TF1* fitFunc = new TF1("fitFunc","crystalball(0) + landau(3)",min0,max0);
    fitFunc->SetParameters(max,maxBinCenter,std,1500,2);
    */

    int success = hist->Fit(fitFunc,"R");

    if (success != 0) {
        cout << "Error: Fit failed with status: " << std::to_string(success) << endl;
    } else{
        cout << "Fit successful" << endl;
    }



    // // Perform the fit
    // TF1* fit = performFit(hist, &max, &maxBinCenter, &std, min0, max0);

    // Drawing the histogram
    TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
    hist->Draw();
    fitFunc->Draw("same");

    // Save and close the file
    c->SaveAs(("./fitData/complexFit/glau/complexFitPeaks_rbf" + std::to_string(rebinFactor) + "_glau.png").c_str());

    delete fitFunc;
    delete c;

    file->Close();
    delete file;
}