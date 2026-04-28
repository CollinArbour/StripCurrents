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
from helpers.accChrg import format_date,loadAccChrgDict
from load_config import config

hvPattern = re.compile(r"hv(\d+)")

# Arguments ------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-d','--date', type=str, help='Measurement date (YYYYMMDD) to process')
    parser.add_argument('--hv',help='Plot all dark rates for this hv value (\'-hv 3600\')')
    parser.add_argument('--2D',help='Plot Dark rate scans vs accumulated charge as heat map')
    parser.add_argument('--strips',help='Plot normalized dark rate noise per strip')
    parser.add_argument('--WG',help='Plot normalized dark rate noise per WG')

    parser.add_argument('--maxRate',type=float,help='Fix the y-axis\' max rate (for comparison)')
    parser.add_argument('--maxHV',type=float,help='Fix the x-axis\' max HV (for comparison)')

    return parser


# Functions ------------------------------------------------
def mkDarkRateScan(base):
    '''
    Collects dark rates at various HV values, and makes summary text

    Args: 
        base: Path to directory containing all hv points dark rate info. This directory contains sub directories for each HV point
    
    Returns:
        Nothing, but writes output to "darkRateScan.txt" file in {base} directory
    '''
    print(f'Doing dark rate scan for:\t{base}')

    # Instantiate out file and write column headers
    outfl = open(f'{base}/darkRateScan.txt','w')

    outfl.write(f'{'HV(V)':<8}')
    for rate in config.ratesOI:
        outfl.write(f'{rate:>8}')
    outfl.write('\n')
    
    # For each HV point, collect their TMB dumps
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

def plotDarkRateScan(base,date,maxRate=None,maxHV=None):
    '''
    Plots dark rates vs HV for a given date of measurement
    
    Args:
        base: Directory of dark rate data, containing "darkRateScan.txt" created by mkDarkRateScan function
        date: Date of measurement, to be used as a title for the plot
    
    Returns:
        Nothing directly, but saves plot to "darkRateScan.png" in {base} directory
    '''
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

    if maxRate:
        ax.set_ylim(top=maxRate)
        ax.set_ylim(0,maxRate)
    if maxHV:
        ax.autoscale(enable=False,axis='x')
        left,_ = ax.get_xlim()
        ax.set_xlim(900,maxHV)
    ax.set_xlabel('HV (V))')
    ax.set_ylabel('Rates (Hz)')
    ax.set_title(f'{date}')
    ax.legend()
    fig.savefig(f'{base}/darkRateScan.png',format='png',dpi=300)
    
def matchRatesToAccChrg(hv,accChrg):
    '''
    Extracts dark rate values for rates specified in config, and matches them to accChrg

    Args:
        hv: HV (V) overwhich scan will be performed
        accChrg: Dictionary of measurement days ("YYYY-MM-DD"), and accChrg value on that day
    
    Returns:
        charge: Array of all accChrg values present
        allRates: 2D array. Each array in the array are rate values of particular rate (e.g. 0ALCT) corresponding to accChrg in charge array
    '''
    charge = []
    allRates = [[] for _ in config.ratesOI]

    for day in config.measurementDays:
        fmtd_date = format_date(day.split('_')[0])
        print(f'Processing: {fmtd_date}')

        fl = glob.glob(f'data/GasGain/data_raw/{day}/data/dark/hv{hv}/*.txt')
        # print(f'\t{fl}')
        
        if len(fl) == 0:
            print(f'\tThere is no run with hv{hv}')
            continue
        elif len(fl) > 1:
            print(f'\tWARNING, more than 1 hv{hv} file (using the first one): \n\t\t{fl}')
        
        rates = hp.parseRates(fl[0],config.ratesOI)

        charge.append(accChrg[fmtd_date])
        for i,rate in enumerate(config.ratesOI):
            allRates[i].append(rates[rate])
        
    return charge,allRates

def plotRateOverTime(hv,charge,allRates):
    for i,rate in enumerate(allRates):
        plt.plot(charge,rate,marker='.',label=config.ratesOI[i])
        plt.axvline(140,color='black',linestyle='--',alpha=0.5)
        plt.axvline(650,color='blue',linestyle='--',alpha=0.15)
        plt.xlabel('Accumulated Charge (mC/cm)')
        plt.ylabel('Rate (Hz)')
        plt.title(f'Dark Rates at HV={hv} V')
        # plt.legend()
        # plt.savefig(f'plots/darkRates/hv{hv}_darkRates_vs_accChrg.png',format='png',dpi=300)

    # for i,malter in enumerate(config.MalterChecks):
    #     plt.axvline(malter[1], color='blue', linestyle='--', alpha=0.1)
    #     # Add label slightly above the plot line (adjust y as needed)
    #     plt.text(malter[1], plt.ylim()[1]*0.95-0.05*i, f"{malter[0].split('_')[0]}", rotation=90,
    #             va='top', ha='right', fontsize=9)
    plt.legend()
    plt.savefig(f'plots/darkRates/hv{hv}_darkRates_vs_accChrg.png',format='png',dpi=300)


# Processing -----------------------------------------------

# Want to watch dark rates grow over HV for each day
    # make text file in data_processed/dark/ DONE
    # plot text file in data_processed/dark/

def main():
    args = build_parser().parse_args()

    # Produce single day's dark rate scan
    if args.date:
        fmtd_date = format_date(args.date)
        base = f'data/GasGain/data_raw/{args.date}_dyn_recupCF4/data/dark/'

        mkDarkRateScan(base)
        plotDarkRateScan(base,fmtd_date,maxRate=args.maxRate,maxHV=args.maxHV)
    
    # Produce dark rates vs acc. charge plot for specifc HV value
    if args.hv:
        print(f'Plotting dark rates at {args.hv} V over time')
        accChrg = loadAccChrgDict()
        charge, allRates = matchRatesToAccChrg(args.hv,accChrg)
        plotRateOverTime(args.hv,charge,allRates)
    
    if args.strips:
        return 0


if __name__ == "__main__":
    main()