#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

import helpers.StripCurrent as hp
import helpers.DataFile as df
import helpers.DataRun as dr
from helpers.functions import quadSum

from load_config import config

# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

### Strategy
# Need to combine with HV scans
#   Dynamically build argument parser (is this the right path? -M --mode option? )
#   This is becuae filtering criteria changes between two cases

# Files to process
parser.add_argument('-f','--file', type=str, help='Name of I vs t signal text file (Kiethly PicoAmmeter format), no paths or extensions')
parser.add_argument('-r','--reference', type=str, help='Reference text file of I vs t (Kiethly PicoAmmeter format)  e.g. 0V or Dark currents')
parser.add_argument('-d','--dataDir', type=str, default='data/StripScans', help='Directory with data, if other than \"data/StripScans\" (for rapid testing)')

# Configuration
parser.add_argument('--hv', type=int, help='HV point to analyze (Strip scan)')
parser.add_argument('-s','--strip', type=int, help='Strip to analyze (HV scan)')
parser.add_argument('-c','--criteria',type=dict, help='Dictionary defining criteria for runs to process e.g.: {"hv":3600,"strips":6}')
parser.add_argument('--source', type=str, help='Source used (e.g. 127Cd, or 90Sr)')

# Action flags
parser.add_argument('--raw', action='store_true', help='Do not process DataRuns, just use datafile as is')
parser.add_argument('--plotRuns', action='store_true', help='Plot each indvidual run as well')
# parser.add_argument('--noFilter',action='store_true',help='Use all runs, no filtering out any runs')
parser.add_argument('--splitRunFile', type=str ,help='directory path to create and write each run of data file into their own file')
parser.add_argument('-v','--verbose', action='store_true', help='Verbose std output')

parser.add_argument('--stripScan', action='store_true', help='Given an HV point, plot a strip current scan')

# Output configuration
parser.add_argument('--outDir', type=str, default='plots/StripScans', help='Directory to save outputs too')
parser.add_argument('--tag', type=str, default='', help='Tag, to append to output file name')


# Functions -----------------------------------------------
def openDataFiles(args):
    # Creating data files objects
    dfs = {
        'signal': None,
        'reference': None
    }

    if args.file:
        dflnm = f'{args.dataDir}/{args.file}.txt'
        print(f'Loading signal file:\n\t{dflnm}')
        dataObject = df.DataFile(args.file)
        dataObject.parseDataFileText(dflnm)
        dfs['signal'] = dataObject
    else:
        print(f'Error, must have signal file!')

    if args.reference:
        rflnm = f'{args.dataDir}/{args.reference}.txt'
        print(f'Loading reference file:\n\t{rflnm}')
        refObject = df.DataFile(args.reference)
        refObject.parseDataFileText(rflnm)
        dfs['reference'] = refObject
    
    return dfs

def plotStripScan(args,dataObjects):
    stripScans = []
    for _,datObj in dataObjects.items():
        if datObj:
            if args.source:
                stripScans.append(datObj.getStripScan())     
            else:
                stripScans.append(datObj.getStripScan(src=False))     

    strips = stripScans[0][0]
    if len(stripScans) == 2:
        if all(stripScans[0][0] == stripScans[1][0]):
            current = stripScans[1][1] - stripScans[0][1]
            error = quadSum(stripScans[1][2],stripScans[0][2])
            plotScan = [stripScans[0][0],cor_I,cor_sig]
        else:
            print('Error! not the same number of strips')
    else:
        current = stripScans[0][1]
        error = stripScans[0][2]
        
    plt.errorbar(strips,current,yerr=error, marker='.',label='Corrected currents')
    plt.xlabel('Strip Number')
    plt.ylabel('Current (nA)')
    plt.legend()
    plt.savefig(f'{args.outDir}/{args.file}_StripScan.png',format='png',dpi=250)


# Running -----------------------------------------------
def main():
    args = parser.parse_args()

    # Creat DataFile objects
    dataObjects = openDataFiles(args)

    ## Refine configurations 
    sort = None
    if args.criteria:
        # Here I would like implement an option to just give criteria that lets you just pick out a specific HV/Strip combination
        # Since this doesn't need args.hv or args.strip I think this needs to be implemented on class side
        # no sorting needed, just filtering (so actually ok at moment, needs class to handle multiple criteria (think it can actually))
        criteria = args.criteria
    elif args.hv or args.strip:
        criteria = {}
        if args.hv:
            # if Doing a strip scan, you need to pick an HV point and sort across strips (to plot from strips 1->12)
            criteria['hv'] = args.hv
            sort = 'strip' 
        elif args.strip:
            # if Doing an HV scan, you need to pick a strip and sort across HV (so you plot as HV rises)
            criteria['strip'] = args.strip
            sort = 'hv'
    else:
        print('Error, you need some sort of criteria')
    
    # Process data objects
    if not args.raw: 
        for run,datObj in dataObjects.items():
            if datObj is None:
                continue

            if criteria:
                datObj.filterRuns(criteria)
            if sort:
                datObj.sortDataRuns(sort)

            dataObjects[run] = datObj
    
    # Processing also needs to include background corrections
        # Need to make sure that df.matching() is dynamic enough to handle matching hv or strips (use sort value?)
    

    ### ACTIONS     (modes if you will)
    # Plot runs if wanted:
    if args.plotRuns:
        outDir = f'{args.outDir}/{args.file}'
        print(f'Writing out to:\n\t{outDir}')
        os.makedirs(outDir,exist_ok=True)

        for _,datObj in dataObjects.items():
            if datObj is None:
                continue
            datObj.describe()
            allRuns = datObj.getDataRuns()
            for run in allRuns:
                run.plotDataRun(save=outDir,tag=args.tag)
    

    if args.stripScan:
        print(f'Performing strip scan:\n\tHV: {args.hv}')
        if args.reference: print(f'\t~ Background corrected ~')
        plotStripScan(args,dataObjects)

    
            

    # # Interpreted arguments (iargs)
    # iargs = process_args(args)



    # # Plotting individual runs:
    # if args.plotRuns:
    #     odir = f'plots/StripScans/{args.outDir}'
    #     os.makedirs(odir,exist_ok=True)
    #     for mdf in mdf_objects:
    #         mdf.describe()

    #         mruns = mdf.getDataRuns()
    #         for run in mruns:
    #             run.plotDataRun(save=odir)

if __name__ == '__main__':
    main()