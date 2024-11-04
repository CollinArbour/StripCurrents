#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

import src.DataFile as df
import src.DataRun as dr
import src.helpers as hp

def matching(src,drk):
    src = np.array(src)
    drk = np.array(drk)
    mtchs,src_idxs,drk_idxs = np.intersect1d(src[0],drk[0],return_indices=True)
    print('\t\tMatched')
    return src[:,src_idxs],drk[:,drk_idxs]

def quadSum(a,b):
    return np.sqrt(a**2 + b**2)

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

run_nms = []

strip_numbers = [2,6]
#strip_numbers = [2,6,10]
ref_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','241002_refMeasures_S10']
plt_nms = ['241014_S2_Plateau','241014_S6_Plateau','241015_S10_Plateau']
pst_nms = ['241029_S2_HVScan03','241029_S6_HVScan03']
additional_ref = ['240820_refMeasures_S7']

mruns = [ref_nms,plt_nms,pst_nms]
mmarks = ['1','2','x']
mcolors = ['blue','green','black']
lims = (10,600)

# Handling Strip Scans
for i,strip in enumerate(strip_numbers):
    criteria = {'strip':strip} 

    print(f'Loading Data and Dark Currents for strip: {strip}')

    for j,run in enumerate(mruns):
        # Load HV Scan with source
        print(f'\tFor run: {run[i]}')
        sdf = df.DataFile(run[i])
        sdf.parseDataFileText(f'./data/HV_Scans/{run[i]}.txt')
        sdf.filterRuns(criteria)
        sdf.sortDataRuns('hv')

        # Load Dark Current runs
        print(f'\tFor dark run: {run[i]}')
        ddf = df.DataFile(f'{run[i]}_dark')
        ddf.parseDataFileText(f'./data/HV_Scans/{run[i]}_dark.txt')
        ddf.filterRuns(criteria)
        ddf.sortDataRuns('hv')

        # Match the data scans to the dark scans
        print('\t\tMatching the HV points between the data and dark scans')

        # Ref Measurements
        src_hvscan = sdf.getHVScan()
        drk_hvscan = ddf.getHVScan(src=False)
        msrc_hvscan,mdrk_hvscan = matching(src_hvscan,drk_hvscan)

        print('\t\tSubtracting background')
        corrected_hvscan = msrc_hvscan[1] - mdrk_hvscan[1]
        #if any(corrected_hvscan < 0):
        #    print('\t\tError More noise than signal!')

        hvscan = [msrc_hvscan[0], corrected_hvscan, quadSum(msrc_hvscan[2],mdrk_hvscan[2])]

        smask = (src_hvscan[0] > lims[0]) * (src_hvscan[0] < lims[1])
        dmask = (drk_hvscan[0] > lims[0]) * (drk_hvscan[0] < lims[1])

    
    #print(pdrkerr)    

    #plt.errorbar(rhv,rI,yerr=np.abs(rerr),linestyle='',label='HV Scan 1 Data',marker='1',color='blue')

    #plt.errorbar(phv,pI,yerr=np.abs(perr),linestyle='',label='HV Scan 2 Data',marker='+',color='green')

    #this creates a differnt set of graphing instructions for lower values
    '''if lims[0] == 0:


        plt.errorbar(rdrkhv,rdrkI,yerr=np.abs(rdrkerr),linestyle='',label='HV Scan 1 Dark',marker='2',color='blue')
        plt.errorbar(pdrkhv,pdrkI,yerr=np.abs(pdrkerr),linestyle='',label='HV Scan 2 Dark',marker='x',color='green')


        plt.yscale('symlog')
        avg_rI = np.mean(rI)
        avg_rdrkI = np.mean(rdrkI)
        avg_pI = np.mean(pI)
        avg_pdrkI = np.mean(pdrkI)

        plt.axhline(y=avg_rI, color='blue', linestyle='--', alpha=0.2)
        plt.axhline(y=avg_rdrkI, color='blue', linestyle='-.', alpha=0.2)
        plt.axhline(y=avg_pI, color='green', linestyle='--', alpha=0.2)
        plt.axhline(y=avg_pdrkI, color='green', linestyle='-.', alpha=0.2)

        plt.text(x=-125, y=avg_rI, s=f'{avg_rI:.2e}', va='center', color='blue')
        plt.text(x=-125, y=avg_rdrkI, s=f'{avg_rdrkI:.2e}', va='center', color='blue')
        plt.text(x=-125, y=avg_pI, s=f'{avg_pI:.2e}', va='center', color='green')
        plt.text(x=-125, y=avg_pdrkI, s=f'{avg_pdrkI:.2e}', va='center', color='green')
    '''
    #plt.plot(rhv,rI,linestyle='-',label='HV Scan 1 Data',marker='',color='blue')
    #plt.plot(phv,pI,linestyle='-',label='HV Scan 2 Data',marker='.',color='green')


    #plt.xlim((0,600))
    #plt.xlim((3400,3800))
    plt.errorbar(src_hvscan[0][smask],src_hvscan[1][smask],yerr=src_hvscan[2][smask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1}')
    plt.errorbar(drk_hvscan[0][dmask],drk_hvscan[1][dmask],yerr=drk_hvscan[2][dmask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1} Dark')

    plt.title(f'Strip {strip} Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.yscale('log')
    plt.tight_layout()
    plt.legend()
    plt.grid(linestyle='-', alpha=0.75, which='both')

    #plt.savefig(f'./plots/HV_Scans/all/full_range/S{strip}_GasGain_HV_Scan2_full_range.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/full_range/S{strip}_GasGain_HV_Scan_full_range.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_GasGain_RefAndPlat_uncorrected_highHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_GasGain_RefAndPlat_uncorrected_plateau_log.png',format='png',dpi=400)

    plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_lowHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_highHV_log.png',format='png',dpi=400)
    plt.close()
