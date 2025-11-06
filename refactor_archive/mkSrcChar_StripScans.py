#!/usr/bin/env python3

from tkinter import font
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

import src.DataFile as df
import src.DataRun as dr
import src.helpers as hp

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)


for i  in range(3):
    if i == 1:
        continue

    fig, ax = plt.subplots(figsize=(8,6))

    flnm = f'./data/StripScans/Src0{i+1}.txt'

    criteria = {'hv':3600}
    
    mdf = df.DataFile(f'Sr_Src0{i+1}')
    mdf.parseDataFileText(flnm)
    mdf.filterRuns(criteria)
    mdf.sortDataRuns('strip')

    mdata = mdf.getStripScan(src=True)
    
    # Convert strip numbers into x position
    x_pos = [0.5*stripWdth]
    for j in range(1,len(mdata[0])):
        x_pos.append(x_pos[-1]+(stripWdth+stripGap))

    # Convert the strip current measurements into linear current densities
    lambda_I = np.array(mdata[1]) / stripWdth   # now has Units now nA/mm
    lambda_I_err = np.array(mdata[2]) / stripWdth   # errors also now nA/mm

    # Format of this data (from DataFile.GetStripScan())
        # mdata[0] = strip numbers
        # mdata[1] = currents
        # mdata[2] = errors

    # strips = [x_pos,mdata[1],mdata[2]]
    strips = [x_pos,lambda_I,lambda_I_err]

    # Fitting the curve to sum of two Gaussians  
    # Initializing parameters with rough estimate
    p0 = [12,11,77,6,24,77] 
    p1, cov = curve_fit(hp.mGaussianSum,np.array(x_pos),strips[1],p0)

    # Convert cross section sigma(x,y=0)
    sig0 = p1[1]
    sig1 = p1[4]
    a0 = p1[0] / (np.sqrt(2*np.pi) * sig0) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = p1[3] / (np.sqrt(2*np.pi) * sig1) # Amplitude of second 2D Gaussian (nA/mm^2)
    px = [a0,sig0,p1[2],a1,sig1,p1[5]]

    # #-------------------------------------------------- 
    # # Make plot showing 2D distribution
    # hp.mkHeatMap_GaussSum(40,p1,mlabel=f'TSrc {i+1}',logScale=False,txtBox=True,save=True)
    # #-------------------------------------------------- 



    #-------------------------------------------------- 
    # Making 1D plots

    # Make plot displaying Strip scan shape and fit
    # ax = hp.mkScans(strips,p1,ax,f'Src{i+1:02}')
    ax = hp.mkScans(strips,px,ax,f'Src{i+1:02}',xsec=True)

    title_fsize = 20
    cntxt = '904 Lab'
    if i==0:
        src ='$^{90}$Sr (30 MBq)' 
        # src ='$^{109}$Cd' 
        # plt.title('$^{90}$Sr (25 MBq) Source',fontsize=title_fsize)
        # plt.title('Linear Current Density for 25 MBq $^{90}$Sr Source')
    elif i==2:
        src ='$^{90}$Sr (50 MBq)' 
        # plt.title('$^{90}$Sr (50 MBq) Source',fontsize=title_fsize)

    lbl_fsize=18
    hp.add_cms_label(plt.gcf(),ax=plt.gca(),cms_fontsize=lbl_fsize,status_fontsize=lbl_fsize,context_fontsize=lbl_fsize,context=cntxt)
    ax.text(0.96,0.96,src,transform=ax.transAxes,fontsize=12,ha='right',va='top')
    plt.gca().tick_params(direction='in')
    plt.legend()

    # plt.savefig(f'./plots/SrSrcs/TSrc{i+1}_avgI.png',format='png',dpi=400)
    # plt.savefig(f'./plots/SrSrcs/TSrc{i+1}_avgI.pdf',format='pdf')

    plt.savefig(f'./plots/SrSrcs/TSrc{i+1}_xSec.png',format='png',dpi=400)
    plt.savefig(f'./plots/SrSrcs/TSrc{i+1}_xSec.pdf',format='pdf')

    plt.close()
    #-------------------------------------------------- 

# NEED TO EDIT THE TITLE!
# plt.title(f'MiniCSC4: Strip Scan L1 HV3400')
# plt.legend()
# plt.savefig(f'./plots/SrSrcs/allSrcs.png',format='png',dpi=400)
# plt.close()
