#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

import helpers.StripCurrent as hp
import helpers.DataFile as df
import helpers.DataRun as dr

from load_config import config

# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

# parser.add_argument('--0VFile', type=str, help='Strip current file at 0V')
parser.add_argument('--darkFile', type=str, help='Dark current txt file (no extension)')
# parser.add_argument('--hv', type=int, help='HV of run to analyze (criteria)')
parser.add_argument('-r','--plotRuns', action='store_true', help='Plot each indvidual run as well')

args = parser.parse_args()

if args.darkFile:
    drkflnm = f'data/StripScans/{args.darkFile}.txt'
    print('Dark file: {drkflnm}')
else:
    print('No Dark File Provided')


# Processing -----------------------------------------------
# fig, ax = plt.subplots(figsize=(8,6))

mdf = df.DataFile(args.darkFile)
mdf.parseDataFileText(drkflnm)
mdf.filterRuns(config.stripCriteria)
mdf.sortDataRuns('strip')

# Plotting individual runs:
if args.plotRuns:
    mdf.describe()
    odir = f'plots/StripScans/{args.darkFile}'
    os.makedirs(odir,exist_ok=True)

    mruns = mdf.getDataRuns()
    for run in mruns:
        run.plotDataRun(save=odir)

# mdata = mdf.getStripScan(src=False)

# ax = hp.mkRawScans(mdata,ax)

# plt.title(f'MiniCSC4: Dark Strip Scan L1 HV3600')
# plt.legend()
# plt.savefig(f'./plots/StripScans/{args.darkFile}.png',format='png',dpi=300)
# plt.close()