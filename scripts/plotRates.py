#!/usr/bin/env python3
import re
import glob
import os
import sys
from turtle import right
import numpy as np
import matplotlib.pyplot as plt
import argparse

import helpers.rates as hp
from helpers.accChrg import format_date
from load_config import config

hvPattern = re.compile(r"hv(\d+)")

# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-d','--date', type=str, help='Measurement date (YYYYMMDD) to process')
parser.add_argument('--dark',action='store_true',help='Do dark rate scan')
parser.add_argument('--hv',help='Plot all dark rates for this hv value (\'-hv 3600\')')
parser.add_argument('--maxRate',type=float,help='Fix the y-axis\' max rate (for comparison)')
parser.add_argument('--maxHV',type=float,help='Fix the x-axis\' max HV (for comparison)')

args = parser.parse_args()

# Functions ------------------------------------------------
def mkDarkRateScan(base):
    print(f'Doing dark rate scan for:\t{base}')

    # Instantiate out file and write column headers
    outfl = open(f'{base}/darkRateScan.txt','w')

    outfl.write(f'{'HV(V)':<8}')
    for rate in config.ratesOI:
        outfl.write(f'{rate:>8}')
    outfl.write('\n')
    
    
    # Process all existing dark rate scans?
    paths = glob.glob(f'{base}/**/*.txt')
    fls = sorted(paths,key=lambda p: os.path.basename(p).lower())

    for fl in fls:
        print(f'\tProcessing: {fl.strip().split('/')[-1]}')
        if match := re.search(hvPattern,fl.lower()):
            hv = match[1]
            outfl.write(f'{hv:<8}')
            
            rates = hp.parseRates(fl,config.ratesOI)
            # Print out collected rates:
            for rate in config.ratesOI:
                outfl.write(f'{rates[rate]:>8.1f}')
            outfl.write('\n')
        else:
            print(f'ERROR: Don\'t recognize file: {fl}')

def plotDarkRateScan(base,date):
    print(f'Making plot for:\t{base}')
    with open(f'{base}/darkRateScan.txt','r') as fl:
        entries = fl.readlines()
    
    lbls = entries[0].strip().split()[1:]
    hv = []
    rates = [[] for _ in lbls]

    for entry in entries[1:]:
        vals = entry.strip().split()
        hv.append(float(vals[0]))

        for i,val in enumerate(vals[1:]):
            rates[i].append(float(val))
    
    fig,ax = plt.subplots()
    
    # print(hv)
    for i,vals in enumerate(rates):
        # print(f'{lbls[i]}: {vals}')
        ax.plot(hv,vals,label=lbls[i],marker='.')

    if args.maxRate:
        ax.set_ylim(top=args.maxRate)
    if args.maxHV:
        ax.autoscale(enable=False,axis='x')
        left,_ = ax.get_xlim()
        ax.set_xlim(3000,args.maxHV)
    ax.set_xlabel('HV (V))')
    ax.set_ylabel('Rates (Hz)')
    ax.set_title(f'{date}')
    ax.legend()
    fig.savefig(f'{base}/darkRateScan.png',format='png',dpi=300)
    


# Processing -----------------------------------------------

# Want to watch dark rates grow over HV for each day
    # make text file in data_processed/dark/ DONE
    # plot text file in data_processed/dark/

if args.dark:
    if not args.date:
        sys.exit("Exit:\tMust include --date argument if doing a dark rate scan")

    base = f'data/GasGain/data_raw/{args.date}_dyn_recupCF4/data/dark/'
    print(f'Doing dark rate scan for:\n\t{base}')
    mkDarkRateScan(base)
    plotDarkRateScan(base,args.date)





# NEED TO REFACTOR THIS
# Want to see dark rates at HV=3600V progression vs accumulated charge
if args.hv:

    accChrg_fl = 'data/GasGain/data_processed/measurement_accChrg.txt'

    # Load Accumulated Charge Values
    with open(accChrg_fl,'r') as fl:
        accChrg_data = fl.readlines()

    accChrg = {}
    for entry in accChrg_data[1:]:
        data = entry.split()
        accChrg[data[0]] = float(data[1])
    
    
    charge = []
    allRates = [[] for _ in config.ratesOI]

    for day in config.measurementDays:
        fmtd_date = format_date(day.split('_')[0])
        print(f'Processing: {fmtd_date}')

        fl = glob.glob(f'data/GasGain/data_raw/{day}/data/dark/hv{args.hv}/*.txt')
        # print(f'\t{fl}')
        
        if len(fl) == 0:
            print(f'\tThere is no run with hv{args.hv}')
            continue
        elif len(fl) > 1:
            print(f'\tWARNING, more than 1 hv{args.hv} file (using the first one): \n\t\t{fl}')
        
        rates = hp.parseRates(fl[0],config.ratesOI)

        charge.append(accChrg[fmtd_date])
        for i,rate in enumerate(config.ratesOI):
            allRates[i].append(rates[rate])
        
    for i,rate in enumerate(allRates):
        plt.plot(charge,rate,marker='.',label=config.ratesOI[i])
        plt.axvline(140,color='black',linestyle='--',alpha=0.5)
        plt.xlabel('Accumulated Charge (mC/cm)')
        plt.ylabel('Rate (Hz)')
        plt.title(f'Dark Rates at HV={args.hv} V')
        plt.legend()
        plt.savefig(f'plots/darkRates/hv{args.hv}_darkRates_vs_accChrg.png',format='png',dpi=300)