#!/usr/bin/env python3

import os
import uproot
import argparse
from scipy.stats import stats,norm
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
#hep.style.use("CMS")

_allDays = ['20240827_dyn_recupCF4',
        '20241025_dyn_recupCF4',
        '20241031_dyn_recupCF4',
        '20241121_dyn_recupCF4',
        '20241122_dyn_recupCF4',
        '20241129_dyn_recupCF4',
        '06_dyn_recupCF4',
        '07_dyn_recupCF4',
        '20250221_dyn_recupCF4'
        ]

parser = argparse.ArgumentParser(description='Fit the charge distribution of a hole to a Gaussian')

parser.add_argument('-d','--days',dest='days',nargs='+',default=_allDays,help='data taking day to fit')
parser.add_argument('-hs','--holes',dest='holes',nargs='*',default=None,help='Name of specific hole to only fit, otherwises does all available holes')
parser.add_argument('-r','--rebin',dest='rebin_factor',type=int,default=64,help='Rebin factor for the charge distribution')
parser.add_argument('-f','--fitWidth',dest='fitWidth',type=float,default=750,help='ChargeADC width around max bin to fit the peak to')
parser.add_argument('-t','--tag',dest='tag',default='',help='Tag to add to the output files')

args = parser.parse_args()

days = args.days

def rebin(hist,rebin_factor):
    if len(hist[0]) % rebin_factor != 0:
        raise ValueError('Number of bins must be divisible by the rebin factor ')
    
    new_counts = np.add.reduceat(hist[0],np.arange(0,len(hist[0]),rebin_factor))
    new_bins = hist[1][::rebin_factor]
    
    return np.histogram(new_bins[:-1],new_bins,weights=new_counts)

def mGauss(x,a,x0,sig):
    return a*np.exp(-(x-x0)**2/(2*sig**2))


# Distribution to fit with Gaussian
dist = 'Cathode/charge/chargeL3'

rebin_factor = args.rebin_factor

# # Original Factors
# rebin_factor = 32
# fitRegion = [1200,2400]

# days = ['07_dyn_recupCF4']



# List of days to fit curves for
for day in days:
    print(f'Fitting {day} runs')

    #################################
    # Get dark Run data and TMB Rates
    dark_path = f'data_processed.nosync/{day}/dark/hv3600'
    darkRunfl = f'{dark_path}/output.root'

    # Open and process dark run
    control = uproot.open(darkRunfl)
    cspec = control[dist].to_numpy()
    rbcspec = rebin(cspec,rebin_factor)
    controlCounts = np.sum(rbcspec[0])

    # Get TMB rate
    with open(f'{dark_path}/TMB_Rate.txt','r') as fl:
        dark_TMBrate = float(fl.readlines()[0].strip().split()[-1])
    print(f'Dark rate for {day} runs is: {dark_TMBrate} Hz')

    #######################################
    # Collect Data Runs and their TMB Rates
    holes = []
    data_path = f'data_processed.nosync/{day}/wSrc'

    wSrcRuns = os.listdir(data_path)
    for run in wSrcRuns:
        if 'hole' in run:
            holes.append(run)

    # Perform fit for each hole
    for hole in holes:
        if args.holes and hole not in args.holes:
            continue
        print(f'\t{hole}')
        dataRunfl = f'{data_path}/{hole}/output.root'
        
        # Open and process data run
        data = uproot.open(dataRunfl)
        dspec = data[dist].to_numpy()
        rbdspec = rebin(dspec,rebin_factor)
        runCounts = np.sum(rbdspec[0])

        # Get TMB rate
        with open(f'{data_path}/{hole}/TMB_Rate.txt','r') as fl:
            run_TMBrate = float(fl.readlines()[0].strip().split()[-1])
        print(f'\t{hole} data run for {day} has TMB rate: {run_TMBrate} Hz')

        # Determine ratio of events due to background in data run
        tmb_ratio = dark_TMBrate/run_TMBrate
        bkg_evts = tmb_ratio*runCounts

        # Get normalized control distribution and make background distribution
        cpdf = rbcspec[0]/np.diff(rbcspec[1]) / controlCounts
        bkg = cpdf * np.diff(rbcspec[1]) * bkg_evts

        # Corrected data distribution
        corrected = (rbdspec[0]-bkg,rbdspec[1])

        # Fit the data--------------------------------------------------
        # Computing bin centers
        centers = 0.5*(corrected[1][1:]+corrected[1][0:-1])

        max_idx = np.argmax(corrected[0])
        maxBin = centers[max_idx]
        fitRegion = [maxBin-args.fitWidth,maxBin+args.fitWidth]

        fitRange = (centers > fitRegion[0]) * (centers < fitRegion[1])
        fitCount = np.sum(corrected[0][fitRange])

        mpdf = corrected[0][fitRange]/32 / fitCount
        mu = np.average(centers[fitRange],weights=mpdf)
        sig = np.sqrt(np.average((centers[fitRange]-mu)**2,weights=mpdf))

        toFit = (mpdf,centers[fitRange])

        p0 = [1,mu,sig]
        p1,cov = curve_fit(mGauss,toFit[1],toFit[0],p0=p0)

        ycheck = mGauss(toFit[1],*p1)
        chi2,pval = stats.chisquare(toFit[0],ycheck)

        xs = np.linspace(fitRegion[0],fitRegion[1],1000)
        ys = mGauss(xs,*p1) * 32 * fitCount

        # Save the fit parameters
        if not args.tag:
            _tag = f'{rebin_factor}_{args.fitWidth}'
        else:
            _tag = args.tag
        
        with open(f'{data_path}/{hole}/fitParams_{_tag}.txt','w') as fl:
            fl.write(f'Peak Position:\t {p1[1]}' + '\n')
            fl.write(f'Sigma:\t\t {p1[2]}' + '\n')
            fl.write(f'Peak Value:\t {p1[0]}' + '\n')
            fl.write(f'Std Error:\t {p1[2]/np.sqrt(fitCount)}' + '\n')
            fl.write(f'Chi Squared: {chi2}'+ '\n')
            fl.write(f'ndof: {len(mpdf)-len(p1)}' + '\n')

        # Plot the graph and save it
        hep.histplot(rbdspec,color='grey',label='Raw Spectrum',alpha=0.8)
        hep.histplot(corrected,color='blue',label='Bkg Subtracted',alpha=0.8)
        plt.plot(xs,ys,color='black',label='Best Fit')

        plt.xlim((0,5000))
        plt.xlabel('Charge in ADC Channels')
        plt.ylabel('Events')
        plt.title(f'109Cd {hole}')

        # Add textbox
        ylims = plt.ylim()
        yval = ylims[0]+0.6*np.diff(ylims)
        plt.text(3000,yval,f'Peak pos: {p1[1]:.1f}\nSigma: {p1[2]:.1f}')

        plt.legend()

        #plt.show()
        plt.savefig(f'{data_path}/{hole}/chargeFit_{_tag}.png',format='png',dpi=200)
        plt.close()