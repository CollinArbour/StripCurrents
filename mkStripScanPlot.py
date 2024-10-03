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
    return src[:,src_idxs],drk[:,drk_idxs]

def quadSum(a,b):
    return np.sqrt(a**2 + b**2)

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','240820_refMeasures_S7','241002_refMeasures_S10']


# Handling Strip Scans
for run_nm in run_nms:
    criteria = {'strip':int(run_nm.split('S')[-1])} 

    # Load HV Scan with source
    sdf = df.DataFile(f'{run_nm[7:]}')
    sdf.parseDataFileText(f'./data/HV_Scans/{run_nm}.txt')
    sdf.filterRuns(criteria)
    sdf.sortDataRuns('hv')

    # Load Dark Current runs
    ddf = df.DataFile(f'{run_nm[7:]}_dark')
    ddf.parseDataFileText(f'./data/HV_Scans/{run_nm}_dark.txt')
    ddf.filterRuns(criteria)
    ddf.sortDataRuns('hv')

    # Match the data points from the two scans
    src_hvscan = sdf.getHVScan()
    drk_hvscan = ddf.getHVScan(src=False)
    src_hvscan,drk_hvscan = matching(src_hvscan,drk_hvscan)

    mhvscan = [src_hvscan[0], src_hvscan[1] - drk_hvscan[1], quadSum(src_hvscan[2],drk_hvscan[2])]

    plt.errorbar(mhvscan[0],mhvscan[1],yerr=mhvscan[2],marker='.',linestyle='')
    plt.savefig('./Source2Test.png',format='png',dpi=400)
    plt.close()


    break
