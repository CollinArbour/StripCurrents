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

run_nms = ['241002_refMeasures_S6']
#run_nms = ['241001_refMeasures_S2','240820_refMeasures_S7','241002_refMeasures_S10']
#run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','240820_refMeasures_S7','241002_refMeasures_S10']


# Handling Strip Scans
for run_nm in run_nms:
    criteria = {'strip':int(run_nm.split('S')[-1])} 
    strip = run_nm.split("_")[-1]

    # Load HV Scan with source
    print(f'Loading Data for {run_nm}')
    sdf = df.DataFile(f'{run_nm[7:]}')
    sdf.parseDataFileText(f'./data/HV_Scans/{run_nm}.txt')
    sdf.filterRuns(criteria)
    sdf.sortDataRuns('hv')

    # Load Dark Current runs
    print(f'Loading Dark Currents for {run_nm}')
    ddf = df.DataFile(f'{run_nm[7:]}_dark')
    ddf.parseDataFileText(f'./data/HV_Scans/{run_nm}_dark.txt')
    ddf.filterRuns(criteria)
    ddf.sortDataRuns('hv')

    # Match the data points from the two scans
    print('\tMatching the HV points of the scan')
    src_hvscan = sdf.getHVScan()
    drk_hvscan = ddf.getHVScan(src=False)
    src_hvscan,drk_hvscan = matching(src_hvscan,drk_hvscan)

    print('\tSubtracting background')
    corrected_hvscan = src_hvscan[1] - drk_hvscan[1]
    if any(corrected_hvscan < 0):
        print('\t\tError More noise than signal!')

    mhvscan = [src_hvscan[0], corrected_hvscan, quadSum(src_hvscan[2],drk_hvscan[2])]

    # Perform fit
    print('\tPerforming fit')
    strt_fit = np.where(mhvscan[0]==3000)[0][0]
    stop_fit = np.where(mhvscan[0]==3550)[0][0]+1
    p0 = [0.001,0.01]
    p1,cov = curve_fit(hp.mExp, mhvscan[0][strt_fit:stop_fit], corrected_hvscan[strt_fit:stop_fit],p0)

    xs = np.linspace(0,3850,2000)
    ys = hp.mExp(xs,p1[0],p1[1])

    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(mhvscan[0],mhvscan[1],yerr=mhvscan[2],marker='.',linestyle='')
    plt.plot(xs,ys)

    plt.axvline(mhvscan[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    plt.axvline(mhvscan[0][stop_fit],color='grey',alpha=0.3,linestyle='--')

    plt.title(f'{strip} Strip Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')

    plt.savefig(f'./plots/HV_Scans/{strip}_RawFitted.png',format='png',dpi=400)
    plt.close()

    # Create a figure
    fig, (ax0, ax1) = plt.subplots(2,1, gridspec_kw={'height_ratios': [4, 3]},sharex=True)
    fig.subplots_adjust(hspace=0)

    # Make Space Charge evaluation
    xs = np.linspace(3000,3850,1500)
    ys = hp.mExp(xs,p1[0],p1[1])

    print('\tCreating Space Charge evaluation')
    ax0.errorbar(mhvscan[0][strt_fit:],mhvscan[1][strt_fit:],yerr=mhvscan[2][strt_fit:],marker='.',linestyle='')
    ax0.plot(xs,ys)
    ax0.set_yscale('symlog')
    ax0.set_title(f'{strip} Strip Current over HV Scan')
    ax0.set_ylabel('Avg. I (nA)')

    yexpect = hp.mExp(mhvscan[0][strt_fit:],p1[0],p1[1])
    diff = yexpect - corrected_hvscan[strt_fit:]
    rel_diff = diff / corrected_hvscan[strt_fit:]

    ax1.plot(mhvscan[0][strt_fit:],rel_diff)
    ax1.set_xlabel('HV (V)')
    ax1.set_ylabel('(Exp-Data)/Data')

    ax0.axvline(mhvscan[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    ax0.axvline(mhvscan[0][stop_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(mhvscan[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(mhvscan[0][stop_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axhline(0,color='black',alpha=0.5,linestyle=':')

    plt.savefig(f'./plots/HV_Scans/{strip}_SpaceCharge.png',format='png',dpi=400)
    plt.close()
