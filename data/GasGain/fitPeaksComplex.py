#!/usr/bin/env python3
import os
os.environ["PYROOT_DISABLE_AUTO_ATEXIT"] = "1"  # 🔥 STOP ROOT's atexit crashes

import ROOT

def getTMBRate(path):
    with open(path, 'r') as fl:
        for line in fl:
            if 'TMB_Rate' in line:
                return float(line.split()[-1])
    return None

def fitPeaksComplex():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1111)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.TH1.AddDirectory(False)

    root_file = './20250401_dyn_recupCF4/wSrc/hv3600/hole2/output.root'
    srcTMBfl = './20250401_dyn_recupCF4/wSrc/hv3600/hole2/TMB_Rate.txt'
    drkFlnm = './20250401_dyn_recupCF4/dark/hv3600/output.root'
    drkTMBfl = './20250401_dyn_recupCF4/dark/hv3600/TMB_Rate.txt'
    hist_name = "Cathode/charge/chargeL3"
    rebin_factor = 128

    srcFile = ROOT.TFile.Open(root_file, "READ")
    drkFile = ROOT.TFile.Open(drkFlnm, "READ")
    if not srcFile or not drkFile:
        print("❌ Could not open one or both ROOT files")
        return

    srcHist = srcFile.Get(hist_name)
    drkHist = drkFile.Get(hist_name)

    if not srcHist or not drkHist:
        print("❌ Missing histograms")
        return

    srcHist = srcHist.Clone("srcHist_clone")
    drkHist = drkHist.Clone("drkHist_clone")
    srcHist.SetDirectory(0)
    drkHist.SetDirectory(0)
    srcHist.Sumw2()
    drkHist.Sumw2()
    srcHist.Rebin(rebin_factor)
    drkHist.Rebin(rebin_factor)

    srcTMB = getTMBRate(srcTMBfl)
    drkTMB = getTMBRate(drkTMBfl)
    if not srcTMB or not drkTMB:
        print("❌ Could not read TMB rates")
        return

    tmb_ratio = drkTMB / srcTMB
    drkHist.Scale(-1 * (tmb_ratio * srcHist.Integral() / drkHist.Integral()))
    correctedHist = srcHist.Clone("correctedHist")
    correctedHist.Add(drkHist)
    correctedHist.SetDirectory(0)

    max_val = correctedHist.GetMaximum()
    max_bin_center = correctedHist.GetBinCenter(correctedHist.GetMaximumBin())
    std = 100.0
    start, stop = 100, 8000

    fitFunc = ROOT.TF1("fitFunc", "crystalball(0) + gaus(5)", start, stop)
    fitFunc.SetParameters(max_val, max_bin_center, std, 3, 2, 100, 3000, 2000)
    ROOT.SetOwnership(fitFunc, False)

    success = correctedHist.Fit(fitFunc, "R")
    if success != 0:
        print(f"❌ Fit failed with status: {success}")
    else:
        print("✅ Fit successful!")

    fitParams = [fitFunc.GetParameter(i) for i in range(8)]

    c = ROOT.TCanvas("canvas_fit", "Fit", 800, 600)
    correctedHist.Draw("HISTE")
    fitFunc.Draw("same")

    cb_part = ROOT.TF1("cb_part", "crystalball(0)", start, stop)
    cb_part.SetParameters(*fitParams[:5])
    cb_part.SetLineColor(ROOT.kRed)
    cb_part.SetLineStyle(2)
    cb_part.Draw("same")

    gaus_part = ROOT.TF1("gaus_part", "gaus(0)", start, stop)
    gaus_part.SetParameters(*fitParams[5:])
    gaus_part.SetLineColor(ROOT.kBlue)
    gaus_part.SetLineStyle(2)
    gaus_part.Draw("same")

    out_dir = "./fitData/complexFit"
    os.makedirs(out_dir, exist_ok=True)
    filename = f"{out_dir}/pyComplexFitPeaks_rbf{rebin_factor}.png"
    c.SaveAs(filename)
    print("📦 Saved:", filename)

    # Manual cleanup
    srcFile.Close()
    drkFile.Close()
    del c, fitFunc, cb_part, gaus_part

# Run it
fitPeaksComplex()
