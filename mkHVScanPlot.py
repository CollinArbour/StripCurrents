#!/usr/bin/env python3
'''

This File is currently being used to process the HV_Scans data files and graph the data using the functions in the helpers class

'''
from re import T
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from itertools import zip_longest

import src.DataFile as df
import src.DataRun as dr
import src.helpers as hp

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

#NOTE: lists of different run files.comment and uncomment as needed. should be in a .data/HV_Scans folder.
run_nms = ['241001_refMeasures_S2', '241002_refMeasures_S6', '241002_refMeasures_S10'] #, '240820_refMeasures_S7'
# run_nms = ['241014_Plateau_S2', '241014_Plateau_S6', '241015_Plateau_S10']

#list for making gasGain tables
table_data = [[], [], []]     #2d list: [[[mhvscan], [holes], [strips]]

# Handling Strip Scans
for run_nm in run_nms:
    criteria = {'strip':int(run_nm.split('S')[-1])} 
    strip = run_nm.split("_")[-1]

    # Load HV Scan with source
    print(f'Loading Data for {run_nm}')
    sdf = hp.createRun(run_nm, criteria)

    # Load Dark Current runs
    print(f'\nLoading Dark Currents for {run_nm}')
    ddf = hp.createRun(run_nm, criteria, src=False)

    # Match the data points from the two scans
    print('\n\tMatching the HV points of the scan')
    src_hvscan = sdf.getHVScan()
    drk_hvscan = ddf.getHVScan(src=False)
    src_hvscan,drk_hvscan = hp.matching(src_hvscan,drk_hvscan)

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
    mhvscan = [src_hvscan[0], corrected_hvscan, hp.quadSum(src_hvscan[2],drk_hvscan[2])]

    src = sdf.getFileSrc()
    hole = sdf.getHole()

    #Generate Raw Fitted Graph
    #hp.mkRawFittedPlot(mhvscan, strip)

    # ---------------------------------------------------------------------------
    # #Generate Gas Gain
    # # hp.mkGasGain(mhvscan, strip, hole,tag='RefMeasures', start_volt=3250,end_volt=3550)
    # hp.mkGasGain(mhvscan, strip, hole,tag='Plateau', start_volt=3250,end_volt=3550)

    # #Gather data for gas gain table making
    # table_data[0].append(mhvscan)
    # table_data[1].append(hole)
    # table_data[2].append(strip)

    # #make gas gain table
    # hp.mkGasGainTable(table_data[0], table_data[1], table_data[2],label=f'plateau_{strip}')
    # ---------------------------------------------------------------------------



    # #create plateau graph
    # hp.mkPlateauPlot(mhvscan, strip, src, hole, uncorrected_curr=src_hvscan, uncorrected_dark_curr=drk_hvscan, y_upper_lim=0.0035)

    #create space charge plot
    hp.mkSpaceChargePlot(mhvscan, strip, start_volt=3300, end_volt=3500)