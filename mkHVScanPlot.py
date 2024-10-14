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
#run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','241002_refMeasures_S10']
#run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','240820_refMeasures_S7','241002_refMeasures_S10']

#run_nms = ['241014_Plateau_S6']
#run_nms = ['241002_refMeasures_S10']

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
        #TODO: add loop that checks for val that triggered error, then create function to show the details of that run
        print('\t\tError More noise than signal!')
        #print("CORRECTED", corrected_hvscan)
        sdf.describe()
        ddf.describe()
    #Holds HV, corrected Avg Current, combined stderror
    mhvscan = [src_hvscan[0], corrected_hvscan, quadSum(src_hvscan[2],drk_hvscan[2])]

    src = sdf.getFileSrc()

    #Generate Raw Fitted Graph(THIS BOTH GENERATES A GRAPH AND RETURNS A VALUE)
    p1 = hp.mkRawFittedPlot(mhvscan, strip)

    #create plateau graph
    sdf.describe()
    ddf.describe()
    hp.mkPlateauPlot(mhvscan, strip, src, end_point=500, uncorrected_curr=src_hvscan, uncorrected_dark_curr=drk_hvscan, y_upper_lim=0.003)

    #create space charge plot
    hp.mkSpaceChargePlot(mhvscan, strip, p1)
