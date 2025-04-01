#!/usr/bin/env python3
import ROOT
import os

def fitPeaksComplex():
    # Enable batch mode (no pop-up windows)
    ROOT.gROOT.SetBatch(True)

    # Enable stats and fit boxes
    ROOT.gStyle.SetOptStat(1111)
    ROOT.gStyle.SetOptFit(1111)

    # File and histogram name
    root_file = "./fitData/output_250221_Hole2_HV3600.root"
    hist_name = "Cathode/charge/chargeL3"

    # Open the ROOT file
    file = ROOT.TFile.Open(root_file, "READ")
    if not file or file.IsZombie():
        print("Error: Could not open file.")
        return

    # Get the histogram
    hist = file.Get(hist_name)
    if not hist:
        print("Error: Could not retrieve histogram.")
        file.Close()
        return

    rebin_factor = 128
    hist.Rebin(rebin_factor)

    # Get initial fit parameters
    max_val = hist.GetMaximum()
    max2 = max_val / 2
    max_bin = hist.GetMaximumBin()
    max_bin_center = hist.GetBinCenter(max_bin)
    std = 100.0
    std2 = 300.0

    fit_range = 800
    min0 = max_bin_center - 1 * fit_range
    max0 = max_bin_center + 2 * fit_range

    # Fit: Gaussian + Landau
    fitFunc = ROOT.TF1("fitFunc", "gaus(0) + landau(3)", min0, max0)
    fitFunc.SetParameters(max2, max_bin_center, std2,
                          max_val, max_bin_center, std)

    # Perform the fit
    success = hist.Fit(fitFunc, "R")

    if success != 0:
        print(f"❌ Fit failed with status: {success}")
    else:
        print("✅ Fit successful!")

    # Draw and save
    c = ROOT.TCanvas("c", "Gaussian Fit", 800, 600)
    hist.Draw()
    fitFunc.Draw("same")

    out_dir = "./fitData/complexFit"
    os.makedirs(out_dir, exist_ok=True)
    filename = f"{out_dir}/pyComplexFitPeaks_rbf{rebin_factor}.png"
    c.SaveAs(filename)

    # Cleanup
    file.Close()

# Run the function
fitPeaksComplex()
