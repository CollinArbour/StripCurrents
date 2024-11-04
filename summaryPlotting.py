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

        plt.errorbar(src_hvscan[0][smask],src_hvscan[1][smask],yerr=src_hvscan[2][smask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1}')
        plt.errorbar(drk_hvscan[0][dmask],drk_hvscan[1][dmask],yerr=drk_hvscan[2][dmask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1} Dark')

    plt.title(f'Strip {strip} Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.yscale('symlog')

    plt.legend()

    plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_lowHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_highHV_log.png',format='png',dpi=400)
    plt.close()
