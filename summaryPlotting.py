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

strip_numbers = [2,6,10]
ref_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','241002_refMeasures_S10']
plt_nms = ['241014_S2_Plateau','241014_S6_Plateau','241015_S10_Plateau']

additional_ref = ['240820_refMeasures_S7']

# Handling Strip Scans
for i,strip in enumerate(strip_numbers):
    criteria = {'strip':strip} 

    # Load HV Scan with source
    print(f'Loading Data for strip: {strip}')
    rsdf = df.DataFile(ref_nms[i])
    psdf = df.DataFile(plt_nms[i])

    rsdf.parseDataFileText(f'./data/HV_Scans/{ref_nms[i]}.txt')
    psdf.parseDataFileText(f'./data/HV_Scans/{plt_nms[i]}.txt')

    rsdf.filterRuns(criteria)
    psdf.filterRuns(criteria)

    rsdf.sortDataRuns('hv')
    psdf.sortDataRuns('hv')

    # Load Dark Current runs
    print(f'Loading Dark Currents for S{strip}')
    rddf = df.DataFile(ref_nms[i])
    pddf = df.DataFile(plt_nms[i])

    rddf.parseDataFileText(f'./data/HV_Scans/{ref_nms[i]}_dark.txt')
    pddf.parseDataFileText(f'./data/HV_Scans/{plt_nms[i]}_dark.txt')

    rddf.filterRuns(criteria)
    pddf.filterRuns(criteria)

    rddf.sortDataRuns('hv')
    pddf.sortDataRuns('hv')

    # Match the data  scans to the dark scans
    print('\tMatching the HV points between the data and dark scans')
    # Ref Measurements
    rsrc_hvscan = rsdf.getHVScan()
    rdrk_hvscan = rddf.getHVScan(src=False)
    rsrc_hvscan,rdrk_hvscan = matching(rsrc_hvscan,rdrk_hvscan)

    # Plateau Measurements
    psrc_hvscan = psdf.getHVScan()
    pdrk_hvscan = pddf.getHVScan(src=False)
    psrc_hvscan,pdrk_hvscan = matching(psrc_hvscan,pdrk_hvscan)

    print('\tSubtracting background')
    rcorrected_hvscan = rsrc_hvscan[1] - rdrk_hvscan[1]
    pcorrected_hvscan = psrc_hvscan[1] - pdrk_hvscan[1]
    #if any(corrected_hvscan < 0):
    #    print('\t\tError More noise than signal!')

    mrhvscan = [rsrc_hvscan[0], rcorrected_hvscan, quadSum(rsrc_hvscan[2],rdrk_hvscan[2])]
    mphvscan = [psrc_hvscan[0], pcorrected_hvscan, quadSum(psrc_hvscan[2],pdrk_hvscan[2])]

    # Plotting all 4 measurements for each strip
    
    # Selecting Data
    lims = (3400,3800)
    rmask = (rsrc_hvscan[0] > lims[0]) * (rsrc_hvscan[0] < lims[1])
    rdrkmask = (rdrk_hvscan[0] > lims[0]) * (rdrk_hvscan[0] < lims[1])

    pmask = (psrc_hvscan[0] > lims[0]) * (psrc_hvscan[0] < lims[1])
    pdrkmask = (pdrk_hvscan[0] > lims[0]) * (pdrk_hvscan[0] < lims[1])

    rhv = rsrc_hvscan[0][rmask]
    rI = rsrc_hvscan[1][rmask]
    rerr = rsrc_hvscan[2][rmask]

    rdrkhv = rdrk_hvscan[0][rdrkmask]
    rdrkI = rdrk_hvscan[1][rdrkmask]
    rdrkerr = rdrk_hvscan[2][rdrkmask]

    phv = psrc_hvscan[0][pmask]
    pI = psrc_hvscan[1][pmask]
    perr = psrc_hvscan[2][pmask]

    pdrkhv = pdrk_hvscan[0][pdrkmask]
    pdrkI = pdrk_hvscan[1][pdrkmask]
    pdrkerr = pdrk_hvscan[2][pdrkmask]
    #print(pdrkerr)    

    plt.errorbar(rhv,rI,yerr=np.abs(rerr),linestyle='',label='HV Scan 1 Data',marker='1',color='blue')
    plt.errorbar(rdrkhv,rdrkI,yerr=np.abs(rdrkerr),linestyle='',label='HV Scan 1 Dark',marker='2',color='blue')

    plt.errorbar(phv,pI,yerr=np.abs(perr),linestyle='',label='HV Scan 2 Data',marker='+',color='green')
    plt.errorbar(pdrkhv,pdrkI,yerr=np.abs(pdrkerr),linestyle='',label='HV Scan 2 Dark',marker='x',color='green')

    #this creates a differnt set of graphing instructions for lower values
    if lims[0] == 0:


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
    
    



    #plt.xlim((0,600))
    #plt.xlim((3400,3800))

    plt.title(f'Strip {strip} Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.yscale('symlog')
    plt.tight_layout()
    plt.legend()

    #plt.show()

    plt.savefig(f'./plots/HV_Scans/all/S{strip}_GasGain_RefAndPlat_uncorrected_highHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_GasGain_RefAndPlat_uncorrected_plateau_log.png',format='png',dpi=400)
    plt.close()
    continue

    # Plotting wholerange (obsolete probably)

    lims = (25,525)

    refplatpts = mrhvscan[1][(mrhvscan[0] > lims[0]) * (mrhvscan[0] < lims[1])]
    pltplatpts = mphvscan[1][(mphvscan[0] > lims[0]) * (mphvscan[0] < lims[1])]

    ref_plat = np.average(refplatpts)
    plt_plat = np.average(pltplatpts)

    print(f'***Strip {strip} ref plat: {ref_plat} \t plateau plat: {plt_plat}')


    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(mrhvscan[0],mrhvscan[1]/ref_plat,yerr=mrhvscan[2]/ref_plat,marker='.',linestyle='',label='Ref Meas.')
    plt.errorbar(mphvscan[0],mphvscan[1]/plt_plat,yerr=mphvscan[2]/plt_plat,marker='.',linestyle='',label='Plt Meas.')

    #plt.plot(xs,ys)
    #plt.axvline(mhvscan[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    #plt.axvline(mhvscan[0][stop_fit],color='grey',alpha=0.3,linestyle='--')

    plt.title(f'{strip} Strip Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.yscale('symlog')

    plt.legend()

    plt.savefig(f'./plots/HV_Scans/all/S{strip}_GasGain_RefAndPlat_log.png',format='png',dpi=400)
    plt.close()
