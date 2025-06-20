void fitSimpleGauss(const char* day, int hv, int hole, int rebinFactor) {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    TString base = "./data_processed/";
    TString hvStr = Form("hv%d", hv);
    TString holeStr = Form("hole%d", hole);
    TString histName = "Cathode/charge/chargeL3";

    TString srcPath = base + day + "/wSrc/" + hvStr + "/" + holeStr;
    TString drkPath = base + day + "/dark/" + hvStr;

    TString srcFileName = srcPath + "/output.root";
    TString drkFileName = drkPath + "/output.root";
    TString srcTMBFile = srcPath + "/TMB_Rate.txt";
    TString drkTMBFile = drkPath + "/TMB_Rate.txt";

    TString outDir = base + day + "/fitPlots";
    gSystem->mkdir(outDir, kTRUE);  // Create directory if it doesn't exist

    TString outPlotName = Form("%s/%s_%s_rbf%d.png", outDir.Data(), hvStr.Data(), holeStr.Data(), rebinFactor);

    // Start opening files
    TFile* srcFile = TFile::Open(srcFileName, "READ");
    TFile* drkFile = TFile::Open(drkFileName, "READ");
    if (!srcFile || srcFile->IsZombie() || !drkFile || drkFile->IsZombie()) {
        printf("❌ Error opening ROOT files\n");
        return;
    }

    // Get histgrams
    TH1D* srcHistOrig = (TH1D*)srcFile->Get(histName);
    TH1D* drkHistOrig = (TH1D*)drkFile->Get(histName);
    if (!srcHistOrig || !drkHistOrig) {
        printf("❌ Error loading histograms\n");
        return;
    }

    // Clone histograms (Isolate them from the files)
    TH1D* srcHist = (TH1D*)srcHistOrig->Clone("srcHistClone");
    TH1D* drkHist = (TH1D*)drkHistOrig->Clone("drkHistClone");
    srcHist->SetDirectory(0);
    drkHist->SetDirectory(0);
    srcFile->Close();
    drkFile->Close();
    delete srcFile;
    delete drkFile;

    srcHist->Rebin(rebinFactor);
    drkHist->Rebin(rebinFactor);

    double drkInt = drkHist->Integral();
    double srcInt = srcHist->Integral();

    auto GetTMBRate = [](const TString& path) {
        ifstream file(path.Data());
        string line;
        while (getline(file, line)) {
            if (line.find("TMB_Rate") != string::npos)
                return stof(line.substr(line.find_last_of(" ") + 1));
        }
        return 1.0f;
    };

    float drkTMB = GetTMBRate(drkTMBFile);
    float srcTMB = GetTMBRate(srcTMBFile);
    float tmb_ratio = drkTMB / srcTMB;

    TH1D* correctedHist = (TH1D*)srcHist->Clone("correctedHist");
    correctedHist->SetDirectory(0);
    correctedHist->Add(drkHist, -1.0 * (1.0 / drkInt) * tmb_ratio * srcInt);

    double maxVal = correctedHist->GetMaximum();
    int maxBin = correctedHist->GetMaximumBin();
    double maxBinCenter = correctedHist->GetBinCenter(maxBin);

    double fitRange = 350;
    double pta = maxBinCenter - fitRange;
    double ptb = maxBinCenter + 1.25*fitRange;

    // Perform Simple Gaussian Fit
    TF1* fitFunc = new TF1("fitFunc", "gaus", pta, ptb);
    fitFunc->SetParameters(maxVal, maxBinCenter, 200);

    int status = correctedHist->Fit(fitFunc, "RQ");
    if (status != 0) printf("❌ Fit failed with status %d\n", status);
    else printf("✅ Fit successful!\n");

    TCanvas* c = new TCanvas("cFit", "Fit Canvas", 800, 600);
    correctedHist->Draw("HIST");
    fitFunc->Draw("same");

    c->SaveAs(outPlotName);
    printf("📦 Output saved: %s\n", outPlotName.Data());

    // Construct output text file name
    TString outTxtName = Form("%s/%s_%s_rbf%d.txt", outDir.Data(), hvStr.Data(), holeStr.Data(), rebinFactor);
    ofstream outFile(outTxtName.Data());

    if (outFile.is_open()) {
        outFile << "Fit Results:\n";
        outFile << "-------------------------\n";
        outFile << Form("Amplitude (p0) = %.4f\n", fitFunc->GetParameter(0));
        outFile << Form("Mean      (p1) = %.4f\n", fitFunc->GetParameter(1));
        outFile << Form("Sigma     (p2) = %.4f\n", fitFunc->GetParameter(2));
        outFile << Form("Chi2             = %.4f\n", fitFunc->GetChisquare());
        outFile << Form("NDF              = %d\n", fitFunc->GetNDF());
        outFile << Form("Chi2/NDF         = %.4f\n", fitFunc->GetChisquare() / fitFunc->GetNDF());
        outFile << Form("Fit Probability  = %.6f\n", fitFunc->GetProb());
        outFile.close();
        printf("📝 Fit summary saved to: %s\n", outTxtName.Data());
    } else {
        printf("❌ Could not open file for writing: %s\n", outTxtName.Data());
    }

    delete c;
    delete fitFunc;
    delete srcHist;
    delete drkHist;
    delete correctedHist;
}
