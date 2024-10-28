#!/usr/bin/env python3
'''

This File is currently being used to process the HV_Scans data files, using the DataFile class, and call upon functions in the helpers class to graph the data.
It also does a small amount of data processing (matching function) and data verifying ("more noise than signal check") before calling upon the graphing functions.

'''
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from itertools import zip_longest

import src.DataFile as df
import src.DataRun as dr
import src.helpers as hp

def matching(src,drk):
    '''
        This function takes in a source scan and dark scan, finds the Matching HV pairs between the two, and removes any nonmatching HV points in either one.

        Arguments:
            -src: hvscan with source
            -drk: hvscan without source
    '''
    src = np.array(src)
    drk = np.array(drk)
    mtchs,src_idxs,drk_idxs = np.intersect1d(src[0],drk[0],return_indices=True)
    print('\t\tMatched')
    return src[:,src_idxs],drk[:,drk_idxs]

def quadSum(a,b):
    return np.sqrt(a**2 + b**2)

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

#NOTE: lists of different run files.comment and uncomment as needed. should be in a .data/HV_Scans folder.
#run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','240820_refMeasures_S7','241002_refMeasures_S10']
#run_nms = ['241014_Plateau_S2', '241014_Plateau_S6', '241015_Plateau_S10']
#run_nms = ['241002_refMeasures_S10']
run_nms = ['241001_refMeasures_S2']

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
    print(f'\nLoading Dark Currents for {run_nm}')
    ddf = df.DataFile(f'{run_nm[7:]}_dark')
    ddf.parseDataFileText(f'./data/HV_Scans/{run_nm}_dark.txt')
    ddf.filterRuns(criteria)
    ddf.sortDataRuns('hv')

    # Match the data points from the two scans
    print('\n\tMatching the HV points of the scan')
    src_hvscan = sdf.getHVScan()
    drk_hvscan = ddf.getHVScan(src=False)
    src_hvscan,drk_hvscan = matching(src_hvscan,drk_hvscan)

    print('\n\tSubtracting background')
    corrected_hvscan = src_hvscan[1] - drk_hvscan[1]

    #Check if the dark measurement is ever larger than the source measurement
    if any(corrected_hvscan < 0):

        #print notfication of error
        print('\n\nERROR: MORE NOISE THAN SIGNAL! PRINTING INDIVIDUAL RUNS!')
        
        #Find the index that triggered the if statement
        index_of_negatives = next((i for i, scan in enumerate(corrected_hvscan) if scan < 0), None)

        #Find the individual run that caused the error to trip, and call describe to print the details of the run
        #NOTE: zip_longest is a imported function that just allows me to loop over both files' dataRuns lists, even if they have different lengths. It does not change anything other than simplify the code.
        for src_run, drk_run in zip_longest(sdf.getDataRuns(), ddf.getDataRuns(), fillvalue='None'):

            if src_run.getHV() == src_hvscan[0][index_of_negatives]:
                sdf.describe(run=src_run)

            if drk_run.getHV() == drk_hvscan[0][index_of_negatives]:
                ddf.describe(run=drk_run)
    
    #Holds HV, corrected Avg Current, combined stderror
    mhvscan = [src_hvscan[0], corrected_hvscan, quadSum(src_hvscan[2],drk_hvscan[2])]

    src = sdf.getFileSrc()
    hole = sdf.getHole()

    #Generate Raw Fitted Graph(THIS BOTH GENERATES A GRAPH AND RETURNS A VALUE [p1])
    p1 = hp.mkRawFittedPlot(mhvscan, strip)

    #Generate Gas Gain
    hp.mkGasGain(mhvscan, strip, start_volt=3000, end_volt=3600)

    #create plateau graph
    #hp.mkPlateauPlot(mhvscan, strip, src, hole, end_point=500, uncorrected_curr=src_hvscan, uncorrected_dark_curr=drk_hvscan, y_upper_lim=0.0035)

    #create space charge plot
    #hp.mkSpaceChargePlot(mhvscan, strip, p1, start_volt=3400, end_volt=3600)
