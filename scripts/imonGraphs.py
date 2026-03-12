#!/usr/bin/env python3

# Imports --------------------------------------------------
import re
import sys
from turtle import end_fill
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from pathlib import Path
import argparse
import pandas as pd 

# import src.helpers as hp
import helpers.accChrg as hp
from load_config import config

irradDays = re.compile(r'(?P<start>\d{6})_(?P<stop>\d{6})')
theDate = re.compile(r"(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})")


# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-d','--date', type=str, default=None, help='Stop date (YYMMDD) of irradiation period')
parser.add_argument('--noAccChrg',action='store_true',help='Don\'t plot accumulated charge vs time')
parser.add_argument('-f','--file', type=str, default=None, help='CAEN log file to plot')
parser.add_argument('-c','--channel', type=int, default=0, help='Which HV channel to read: default = 0 for irradiation')
parser.add_argument('-o','--outfile', type=str, default=None, help='Outfile (path + ) name')
parser.add_argument('-p','--pressure',nargs='?',const='data/Irradiation/climate/pressure.csv',default=None,help='Match and plot pressure data with imon')
parser.add_argument('-t','--temperature',nargs='?',const='data/Irradiation/climate/temperature.csv',default=None,help='Match and plot temperature data with imon')

args = parser.parse_args()

# Functions ------------------------------------------------
def findIrradLog(args):
    # find log file
    logfl = ''
    for fl in config.Period01_src1_log_files + config.Period02_src3_log_files:
        days = re.search(irradDays,fl)
        # print(f'{days.group('start')}, {days.group('stop')}, {args.date}, {(days.group('stop') == args.date)}')
        # print(f'{days.group("nstop")} == {args.date} = (days.group('stop') == args.date)')
        if days and days.group('stop') == args.date:
            assert Path(fl).exists(), f'Missing irradiation log: {fl}'
            logfl = fl

    if not logfl:
        sys.exit(f'ERROR: Configuration is unaware of any irradiation periods ending on: {args.date}')

    dates = re.findall(theDate,logfl)
    start_date = f'20{'-'.join(dates[0])}'
    end_date = f'20{'-'.join(dates[1])}'

    print(start_date, ' to ', end_date)
    
    return logfl,start_date,end_date

def filterTimeRange(df,column,start_time,stop_time):
    df = df.sort_values(column)
    return df.loc[
        (df[column] >= start_time)*(df[column] <= stop_time)
    ].reset_index(drop=True)

def plotCombined(df,sdate,edate,ofl=None):
    lbl_dict = {
        'current':'Current (µA)',
        'temperature':'Temperature (°C)',
        'pressure':'Pressur (mbar)'
    }

    columns = df.columns
    print(columns)
    if len(columns) % 2 != 0:
        print('Odd number of columns provided')
        return

    nPlots = int(len(columns)/2)
    fig, ax = plt.subplots(nPlots,1,sharex=True,squeeze=False)
    ax = ax.flatten()

    for i in range(nPlots):
        j = i*2     # to get to first index of pair of two columns
        ax[i].plot(df[columns[j]],df[columns[j+1]])
        ax[i].set_ylabel(lbl_dict[columns[j+1]])
    
    ax[-1].tick_params(axis='x',labelrotation=25)
    ax[-1].set_xlabel('Time (Date)')
    
    fig.tight_layout()
    
    if ofl:
        plt.savefig(ofl,format='png',dpi=200)
    else:
        plt.savefig(f'plots/accCharge/IMon/{sdate.strftime('%Y%m%d')}_{edate.strftime('%Y%m%d')}_IMon.png', format='png', dpi=200)


# Processing -----------------------------------------------
def main():
    # Identify current file
    if args.date and not args.file:
        print('Processing official irradiation period')
        logfl,start_date,end_date = findIrradLog(args)
    if args.file and not args.date:
        print('Processing general CAEN log file')
        logfl = args.file
        start_date = end_date = None

    # Reads and cleans log data file using regex to only contain timestamps and imon values
    vmatches,imatches = hp.parse_log(logfl,ch=args.channel)

    # Creates list of timestamps and imon values from cleaned log file
    timestamps = hp.timestamps(imatches)
    imon_values = hp.values(imatches)

    idf = pd.DataFrame({
        'Itime': timestamps,
        'current': imon_values
    })

    # Set date range if not provided 
    if start_date is None:
        start_date = idf['Itime'].iloc[0]
    if end_date is None:
        end_date = idf['Itime'].iloc[-1]

    # Filter incase imon logging quit before end of irradiation
    # Need the dates for approximate plotting
    print('Filtering Current')
    idf_filtered = filterTimeRange(idf,'Itime',start_date,end_date)

    idf_resampled = idf_filtered.resample('5min',on='Itime').mean(numeric_only=True).reset_index()

    # Add pressure or temperature data if present
    if args.pressure:
        # NEED TO IMPROVE, DETECT GIF++ OR B904 DATA FORMAT AND TREAT APPROPRIATELY
        print(f'Processing pressure data: {args.pressure}')
        pdf = hp.readPressure(args.pressure)
        pdf_filtered = filterTimeRange(pdf,'Ptime',start_date,end_date)
        idf_resampled = pd.concat([idf_resampled,pdf_filtered],axis=1)
    if args.temperature:
        print(f'Processing temperature data: {args.temperature}')
        tdf = hp.readTemperature(args.temperature)
        tdf_filtered = filterTimeRange(tdf,'Ttime',start_date,end_date)
        idf_resampled = pd.concat([idf_resampled,tdf_filtered],axis=1)

    # # Maybe don't need this? always suspiciously linear
    # if not args.noAccChrg:
    #     print('Also plotting accumulated charge vs time')
    #     #calculates accumulated charge based on imon values and timestamps
    #     accumulated_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')
    #     ''' shows graph of the accumulated charge vs time '''
    #     hp.accCharge_vs_time(start_date, end_date,timestamps, accumulated_charge)

    # ''' shows a graph of the current vs time '''
    # hp.current_vs_time(start_date, end_date, timestamps, imon_values)

    print(idf_resampled)
    plotCombined(idf_resampled,start_date,end_date,ofl=args.outfile)


if __name__ == '__main__':
    main()