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

# parser.add_argument('-fls','--files', type=str, nargs='+', help='Files to process')
parser.add_argument('--noHVfile', type=str, help='Strip current file at 0V (no file extension)')
parser.add_argument('--darkfile', type=str, help='Dark current txt file (no file extension)')
parser.add_argument('--hv', type=int, help='HV of run to analyze (criteria)')
parser.add_argument('-p','--plotRuns', action='store_true', help='Plot each indvidual run as well')
parser.add_argument('--outDir', type=str, required=True, help='Directory to save outputs too')
parser.add_argument('--tag', type=str, default='', help='Tag, to append to output file name')
parser.add_argument('--noFilter',action='store_true',help='Dont\'t reject duplicates, plot all runs')
parser.add_argument('--splitRunFile', action='store_true',help='Print each individual run into their own text file')

args = parser.parse_args()

if args.darkfile and not args.hv:
    parser.error("Need to provide HV of dark runs when --darkfile is used")

toProcess = []
allCriteria = []
stripScans = []

if args.noHVfile:
    noHVflnm = f'data/StripScans/{args.noHVfile}.txt'
    print(f'0V file: {noHVflnm}')
    toProcess.append(noHVflnm)
    allCriteria.append({'hv':0})

if args.darkfile:
    drkflnm = f'data/StripScans/{args.darkfile}.txt'
    print(f'Dark file: {drkflnm}')
    toProcess.append(drkflnm)
    allCriteria.append({'hv':args.hv})
else:
    print('No Files Provided')


# Processing -----------------------------------------------
# fig, ax = plt.subplots(figsize=(8,6))

mdf_objects = []

for i,fl in enumerate(toProcess):
    mdf = df.DataFile(fl)
    mdf.parseDataFileText(fl)
    if not args.noFilter:
        mdf.filterRuns(allCriteria[i])
    mdf.sortDataRuns('strip')
    mdf_objects.append(mdf)

    stripScans.append(mdf.getStripScan(src=False))

# Plotting individual runs:
if args.plotRuns:
    odir = f'plots/StripScans/{args.outDir}'
    os.makedirs(odir,exist_ok=True)
    for mdf in mdf_objects:
        mdf.describe()

        mruns = mdf.getDataRuns()
        for run in mruns:
            run.plotDataRun(save=odir)

# Split up run file into individual runs
if args.splitRunFile:
    for mdf in mdf_objects:
        outDir = f'plots/StripScans/{args.outDir}'
        mdf.splitRuns(outDir,tag=args.tag) 

# Plot the strip scan
if len(stripScans) == 2:
    if all(stripScans[0][0] == stripScans[1][0]):
        cor_I = stripScans[1][1] - stripScans[0][1]
        cor_sig = quadSum(stripScans[1][2],stripScans[0][2])
    
    plt.errorbar(stripScans[0][0],cor_I,yerr=cor_sig, marker='.',label='Corrected currents')
    plt.xlabel('Strip Number')
    plt.ylabel('Current (nA)')
    plt.legend()
    
    plt.savefig(f'plots/StripScans/{args.outDir}/stripScan.png',format='png',dpi=200)