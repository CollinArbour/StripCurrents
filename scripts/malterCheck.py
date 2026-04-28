#!/usr/bin/env python3

# Imports --------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import argparse
from pathlib import Path

import helpers.accChrg as hp
from load_config import config

# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-d','--date', type=str, help='Date (YYMMDD) of Malter current check')
parser.add_argument('-s','--smooth', action='store_true', help='Plot a rolling average of current')
parser.add_argument('-a','--all', action='store_true', help='Plot all Malter checks on the same plot')


# Functions ------------------------------------------------
def processMalter(Vmon,Imon):
    vtime = np.array(hp.timestamps(Vmon))
    itime = np.array(hp.timestamps(Imon))

    if vtime[0] < itime[0]:
        offset = vtime[0]
    else:
        offset = itime[0]

    vtime -= offset
    itime -= offset

    vtime = [dt.total_seconds() for dt in vtime]
    itime = [dt.total_seconds() for dt in itime]

    voltage = pd.Series(data=hp.values(Vmon),index=vtime)
    current = pd.Series(data=hp.values(Imon),index=itime)

    Ismooth = current.rolling(window=10,center=True).mean()

    return voltage, current, Ismooth

def findNoSource(Imon):
    '''Finds the moment in the file where source was removed'''
    return 0


def plotMalter(date,Vmon,Imon,smooth=False):
    voltage, current, Ismooth = processMalter(Vmon,Imon)

    fig, ax = plt.subplots()

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Vmon (V)',color='tab:blue')
    ax.set_ylim(0,3750)
    ax.plot(voltage.index,voltage.values,linestyle='',marker='.',color='tab:blue',label='VMon')
    # ax.plot(voltage.index,voltage.values,color='tab:blue',label='VMon')
    ax.tick_params(axis='x',rotation=25)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(hp.format_mmss))

    ax.tick_params(axis='y',labelcolor='tab:blue')
    ax.ticklabel_format(axis='y',style='plain')

    ax1 = ax.twinx()
    if smooth:
        ax1.plot(Ismooth.index,Ismooth.values,color='tab:red',label='IMon (smooth)')
    else:
        stable = current[current.diff().abs() < 0.03]   # threshold in µA
        ax1.plot(stable.index,stable.values,color='tab:red',label='IMon')
        # ax1.plot(stable.index,stable.values,color='tab:red',linestyle='',marker='.',label='IMon')
        # ax1.plot(current.index,current.values,color='tab:red',linestyle='',marker='.',label='IMon')
        # ax1.plot(current.index,current.values,color='tab:red',label='IMon')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Imon ($\mu A$)', color='tab:red')
    ax1.tick_params(axis='y',labelcolor='tab:red')

    plt.title('Malter Current Check')
    fig.tight_layout()

    outFl = f'./plots/accCharge/Malter/{date}_malter_log'
    if smooth:
        outFl += '_smooth'
    
    plt.savefig(f'{outFl}.png',format='png',dpi=250)
    plt.close()

def plotAllMalter(smooth=False):
    mbase = f'./data/Irradiation/MalterCurrents'

    fig, ax = plt.subplots()

    for check,accChrg in config.MalterChecks:
        MalterFile = Path(f'{mbase}/{check}')
        assert MalterFile.exists(), f'Missing path: {MalterFile}'
        print(f'Processing: {check}')

        vmatches ,imatches = hp.parse_log(MalterFile)
        _ , current, Ismooth = processMalter(vmatches,imatches)

        if smooth:
            to_plot=Ismooth
            outFl = 'allMalterChecks_smooth'
        else:
            to_plot=current
            outFl = 'allMalterChecks'
        
        ax.plot(to_plot.index,to_plot.values,label=f'{accChrg} mC/cm')
    
    ax.tick_params(axis='x',rotation=25)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(hp.format_mmss))
    ax.set_xlabel('Time (min)')

    ax.set_ylabel(r'Imon ($\mu A$)')
    ax.set_yscale('log')

    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f'./plots/accCharge/Malter/{outFl}.png',format='png',dpi=300)


# Processing -----------------------------------------------
if __name__ == '__main__':
    args = parser.parse_args()

    if args.date:
        print('NOT SUPPOSED TO ENTER AT THE MOMENT (developing multi-Malter)')
        MalterFile = Path(f'./data/Irradiation/MalterCurrents/Malter_{args.date}.log')

        # Check that file exists
        assert MalterFile.exists(), f'Missing path: {MalterFile}'

        #reads and cleans log data file using regex to only contain timestamps and imon values
        vmatches,imatches = hp.parse_log(MalterFile)

        print(f'Plotting: {MalterFile}')
        plotMalter(args.date,vmatches,imatches,smooth=args.smooth)

    if args.all:
        print('')
        plotAllMalter(smooth=args.smooth)