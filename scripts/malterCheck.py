#!/usr/bin/env python3

# Imports --------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime
import argparse
from pathlib import Path

import helpers.accChrg as hp
# from load_config import config

# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-d','--date', type=str, required=True, help='Date (YYMMDD) of Malter current check')


# Functions ------------------------------------------------
def plotMalter(date,Vmon,Imon):
    # print('hello world')
    vtime = np.array(hp.timestamps(Vmon))
    vval = hp.values(Vmon)
    
    itime = np.array(hp.timestamps(Imon))
    ival = np.array(hp.values(Imon))


    if vtime[0] < itime[0]:
        offset = vtime[0]
    else:
        offset = itime[0]
    
    vtime -= offset
    itime -= offset

    vtime = [dt.total_seconds() for dt in vtime]
    itime = [dt.total_seconds() for dt in itime]

    fig, ax = plt.subplots()

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Vmon (V)',color='tab:blue')
    ax.plot(vtime,vval,color='tab:blue',label='VMon')
    ax.tick_params(axis='x',rotation=25)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(hp.format_mmss))

    ax.tick_params(axis='y',labelcolor='tab:blue')
    ax.ticklabel_format(axis='y',style='plain')

    ax1 = ax.twinx()
    ax1.plot(itime,ival,color='tab:red',label='IMon')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Imon ($\mu A$)', color='tab:red')
    ax1.tick_params(axis='y',labelcolor='tab:red')

    plt.title('Malter Current Check')
    fig.tight_layout()

    plt.savefig(f'./plots/accCharge/Malter/{date}_malter_log.png',format='png',dpi=250)
    plt.close()



# Processing -----------------------------------------------
if __name__ == '__main__':
    args = parser.parse_args()

    MalterFile = Path(f'./data/Irradiation/MalterCurrents/{args.date}_Malter.log')

    # Check that file exists
    assert MalterFile.exists(), f'Missing path: {MalterFile}'

    #reads and cleans log data file using regex to only contain timestamps and imon values
    vmatches,imatches = hp.parse_log(MalterFile)

    print(f'Plotting: {MalterFile}')
    plotMalter(args.date,vmatches,imatches)
