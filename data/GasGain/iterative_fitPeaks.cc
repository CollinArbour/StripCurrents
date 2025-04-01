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
    gStyle->SetOptStat(0); // Turn off the stats box
    //gStyle->SetOptStat(1111);  // Show entries, mean, and RMS
    //gStyle->SetOptFit(1111);   // Show fit results (Chi2, NDF, parameters)

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

    // Create storage arrays for the fit parameters
    std::vector<double> Ranges = {};
    std::vector<double> Amps = {};
    std::vector<double> Means = {};
    std::vector<double> Sigmas = {};
    std::vector<double> Chi2s = {};
    std::vector<double> NDFs = {};
    std::vector<double> Chi2OverNDF = {};

    // Begin iterative fitting process 
    for (double fitRange = 300; fitRange < 1800; fitRange+=32){
        // Update the fit range
        double min0 = maxBinCenter - fitRange;
        double max0 = maxBinCenter + fitRange;

        cout << "----------------------------------------------------------------------------------------------------" << endl;
        cout << "Peaked at " << maxBinCenter << " with amplitude " << max << endl;
        cout << "Fitting the range from " << min0 << ", to " << max0 <<  " (Range of " << fitRange << ")" << endl;
        cout << endl;
        Ranges.push_back(2*fitRange);

        // Perform the fit
        TF1* fit = performFit(hist, &max, &maxBinCenter, &std, min0, max0);

        // Store the fit parameters
        Amps.push_back(max);
        Means.push_back(maxBinCenter);
        Sigmas.push_back(std);
        Chi2s.push_back(fit->GetChisquare());
        NDFs.push_back(fit->GetNDF());
        Chi2OverNDF.push_back(fit->GetChisquare()/fit->GetNDF());
        
        cout << "GOF: " << fit->GetChisquare()/fit->GetNDF() << endl;

        // Drawing the histogram
        TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
        hist->Draw();
        fit->Draw("same");

        // Save and close the file
        c->SaveAs(("./fitData/rbf128_expandingFitRange/FR.png" + std::to_string(static_cast<int>(fitRange)) + ".png").c_str());

        delete fit;
        delete c;

        cout << maxBinCenter << endl;
        //break;
    }

    // Draw the fit parameters
    TGraph* peakPos_graph = new TGraph(Ranges.size(), Ranges.data(), Means.data());
    // Customize graph
    peakPos_graph->SetTitle("Peak Postiion vs FitRange;Fit Range (ADC);Peak Position");  // Title & axis labels
    peakPos_graph->SetMarkerStyle(20);  // Square markers
    peakPos_graph->SetMarkerSize(1.2);
    peakPos_graph->SetMarkerColor(kBlue);
    // Draw graph
    TCanvas* c1 = new TCanvas("c1", "PeakPos Plot", 800, 600);
    peakPos_graph->Draw("AP");  // "A" (axis), "P" (points)
    // Save plot
    c1->SaveAs("./fitData/peakPosition_vs_fitRange.png");
    delete c1;
    
    // Draw the Chi2/NDF
    TGraph* GOF_graph = new TGraph(Ranges.size(), Ranges.data(), Chi2OverNDF.data());
    GOF_graph->SetTitle("Chi2/NDF vs FitRange;Fit Range (ADC);Chi2/NDF");  // Title & axis labels
    GOF_graph->SetMarkerStyle(20);  // Square markers
    GOF_graph->SetMarkerSize(1.2);
    GOF_graph->SetMarkerColor(kGreen);

    TCanvas* c2 = new TCanvas("c2", "GOF Plot", 800, 600);
    GOF_graph->Draw("AP");  // "A" (axis), "P" (points)
    c2->SaveAs("./fitData/gof_vs_fitRange.png");
    delete c2;
    
    /*
    // Drawing the histogram
    TCanvas* c = new TCanvas("c", "Gaussian Fit", 800, 600);
    hist->Draw();
    mfit->Draw("same");

    // Retrieve the stats box
    gPad->Update();  // Required to make sure the stats box is drawn

    // Save and close the file
    c->SaveAs("./fitData/mfunc_fitPeaks_rbf128_750.png");
    */

    file->Close();

    delete file;
    //delete c;
}