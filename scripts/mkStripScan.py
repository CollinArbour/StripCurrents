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

# Methods to dynamically parse input --------------------------------------------
def add_CommonArgs(parser):
    # Files to process
    parser.add_argument('-f','--file', type=str, help='Name of I vs t signal text file (Kiethly PicoAmmeter format), no paths or extensions')
    parser.add_argument('-r','--reference', type=str, help='Reference text file of I vs t (Kiethly PicoAmmeter format)  e.g. 0V or Dark currents')

    # Action flags
    parser.add_argument('--raw', action='store_true', help='Do not process DataRuns, just use datafile as is')
    # parser.add_argument('--unfiltered', action='store_true', help='Do not filter (sort) DataRuns, just use datafile as is')
    parser.add_argument('--plotRuns', action='store_true', help='Plot each indvidual run as well')
    parser.add_argument('--splitRunFile', type=str ,help='Path to dir, script creates and writes each run of data file into their own file here')
    parser.add_argument('-v','--verbose', action='store_true', help='Verbose std output')

    # Configuration
    parser.add_argument('--source', type=str, help='Source used (e.g. 127Cd, or 90Sr)')

    # Output configuration
    parser.add_argument('--tag', type=str, default='', help='Tag, to append to output file name')

def add_stripScanParser(subparsers):
    parser = subparsers.add_parser(
        'StripScan',
        help = 'Strip scan at fixed HV'
    )

    add_CommonArgs(parser)

    # Configuration
    parser.add_argument('--hv', type=int, default=3600, help='HV point to analyze')
    parser.add_argument('-d','--dataDir', type=str, default='data/StripScans', help='Directory with data, if other than \"data/StripScans\" (for rapid testing)')
    parser.add_argument('--outDir', type=str, default='plots/StripScans', help='Directory to save outputs too (default=plots/StripScans)')
    parser.set_defaults(
        sort = 'strip'
    )

def add_hvScanParser(subparsers):
    parser = subparsers.add_parser(
        'HVScan',
        help = 'HV scan at specific strip'
    )

    add_CommonArgs(parser)

    parser.add_argument('-s','--strip', type=int, help='Strip to analyze (HV scan)')
    parser.add_argument('-d','--dataDir', type=str, default='data/HVScans', help='Directory with data, if other than \"data/HVScans\" (for rapid testing)')
    parser.add_argument('--outDir', type=str, default='plots/HVScans', help='Directory to save outputs too (default=plots/HVScans)')
    parser.set_defaults(
        sort = 'hv'
    )

# Would like to add a mode that just processes general files
#       (general mode, e.g. for just plotting all runs or breaking apart and stitching back together runs) 

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
    # Parse options dynamically
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='mode', required=True)

    add_stripScanParser(subparsers)
    add_hvScanParser(subparsers)

    args = parser.parse_args()


    # Creat DataFile objects -------------------------
    dataObjects = openDataFiles(args)

    ## Refine configurations 
    criteria = {}
    if hasattr(args,'hv'):
        criteria['hv'] = args.hv
    if hasattr(args,'strip'):
        criteria['strip'] = args.strip

    if args.verbose:
        print(f'\tCriteria: {criteria}')
    
    # Process data objects
    if not args.raw: 
        for run,datObj in dataObjects.items():
            if datObj is None:
                print(f'No data object associated with run: {run}')
                continue

            if criteria:
                datObj.filterRuns(criteria)
            if args.sort:
                datObj.sortDataRuns(args.sort)

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
    

    if args.mode == 'StripScan':
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