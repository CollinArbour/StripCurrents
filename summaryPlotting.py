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
lims = (0,3800)
'''
important lims:
    low hv: (0,600)
    high hv: (3400,3800)
    full:(0,3800)

    change above variable to change type
'''


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


        if lims[0] == 0 and lims[1] == 600:
            '''
            this is the low HV scan 
            needs the current values and the dark rate values to be noted on here
            '''
            plt.errorbar(src_hvscan[0][smask],src_hvscan[1][smask],yerr=src_hvscan[2][smask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1}')
            plt.errorbar(drk_hvscan[0][dmask],drk_hvscan[1][dmask],yerr=drk_hvscan[2][dmask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1} Dark')
            
            avg_hvscan = np.mean(src_hvscan[1][smask])
            avg_hvscan_drk = np.mean(drk_hvscan[1][dmask])

            plt.axhline(y=avg_hvscan, color=mcolors[j], linestyle='--', alpha=0.2)
            plt.axhline(y=avg_hvscan_drk, color=mcolors[j], linestyle='-.', alpha=0.2)
            
            plt.text(x=-125, y=avg_hvscan, s=f'{avg_hvscan:.2e}', va='center', color=mcolors[j])
            plt.text(x=-125, y=avg_hvscan_drk, s=f'{avg_hvscan_drk:.2e}', va='center', color=mcolors[j])

            plt.yscale('symlog')
            plt.xlim(left=0)


        if lims[0] == 3400 and lims[1] == 3800:
            '''
            this is the hv scan
            dark rate values are not necessary since they are so small
            '''
            plt.errorbar(hvscan[0],hvscan[1],yerr=hvscan[2],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1}')
            plt.xlim(lims[0]-25,lims[1]+25 )
            plt.grid(linestyle='-', alpha=0.75)


        if lims[0] == 0 and lims[1] == 3800:
            '''
            this is the full ranged scan
            will be plotted in a true log scale with log grid lines, dark rate not needed but use the corrected values
            find a way to use corrected_hvscan
            run matching for corrected hv and src_hvscan[1][smask]
            '''
            #plt.plot(hvscan[0], hvscan[1],linestyle='-', color=mcolors[j], label=f'Scan {j+1}')
            #plt.errorbar(hvscan[0],hvscan[1],yerr=hvscan[2],linestyle='',color=mcolors[j])

            #plt.errorbar(src_hvscan[0][smask],src_hvscan[1][smask],yerr=src_hvscan[2][smask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1}')
            #plt.errorbar(drk_hvscan[0][dmask],drk_hvscan[1][dmask],yerr=drk_hvscan[2][dmask],linestyle='',marker=mmarks[j],color=mcolors[j],label=f'Scan {j+1} Dark')
            
            plt.plot(src_hvscan[0][smask],src_hvscan[1][smask],linestyle='-',color=mcolors[j],label=f'Scan {j+1}')
            plt.plot(drk_hvscan[0][dmask],drk_hvscan[1][dmask],linestyle='-',color=mcolors[j],label=f'Scan {j+1} Dark')
            


            plt.yscale('log')
            plt.grid(linestyle='-', alpha=0.75, which='both')
            plt.xlim(left=0)
            
        
    plt.title(f'Strip {strip} Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.tight_layout()
    plt.legend()
    plt.show()

    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_lowHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_highHV_log.png',format='png',dpi=400)
    #plt.savefig(f'./plots/HV_Scans/all/S{strip}_hvScan_full_range_log.png',format='png',dpi=400)

    