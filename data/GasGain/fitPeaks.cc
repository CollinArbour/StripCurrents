#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>


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
    
    // Find the max value of the histogram
    double max = hist->GetMaximum();
    int maxBin = hist->GetMaximumBin();
    double maxBinCenter = hist->GetBinCenter(maxBin);
    
    // Initial fit to the function
    int fitRange = 300;
    double min0 = maxBinCenter - fitRange;
    double max0 = maxBinCenter + fitRange;
    cout << "--------------------------------------------------" << endl;
    cout << "*** Initial Guess *** " << endl;
    cout << "Max bin center: " << maxBinCenter << endl;
    cout << "Fitting range: " << min0 << " - " << max0 << endl;
    cout << "--------------------------------------------------" << endl;
    TF1* fitFunc1 = new TF1("fitFunc1", "gaus", min0, max0);

    fitFunc1->SetParameters(0,max);
    fitFunc1->SetParameter(1,maxBinCenter);
    fitFunc1->SetParameter(2,fitRange);

    hist->Fit(fitFunc1,"R");

    // Get initial fit parameters
    double chi21 = fitFunc1->GetChisquare();
    int ndf1 = fitFunc1->GetNDF();
    double const1 = fitFunc1->GetParameter(0);
    double mean1 = fitFunc1->GetParameter(1);
    double sigma1 = fitFunc1->GetParameter(2);

    int min1 = mean1 - 1.5*sigma1;
    int max1 = mean1 + 1.5*sigma1;

    cout << "--------------------------------------------------" << endl;
    cout << "*** Initializing 2nd Fit ***" << endl;
    cout << "max: " << const1 << endl;
    cout << "Max bin center: " << mean1 << endl;
    cout << "Fitting range: " << min1 << " - " << max1 << endl;
    cout << "--------------------------------------------------" << endl;

    TF1* fitFunc2 = new TF1("fitFunc2", "gaus", min1, max1);
    fitFunc2->SetParameters(0,max1);   
    fitFunc2->SetParameter(1,mean1);
    fitFunc2->SetParameter(2,sigma1);

    hist->Fit(fitFunc2,"R");

    // Drawing the histogram
    TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
    hist->Draw();
    fitFunc1->Draw("same");
    fitFunc2->Draw("same");

    // Retrieve the stats box
    gPad->Update();  // Required to make sure the stats box is drawn

    // Save and close the file
    c->SaveAs("./fitData/fitPeaks_rbf128.png");
    file->Close();

    delete file;
    delete fitFunc1;
    delete c;
}