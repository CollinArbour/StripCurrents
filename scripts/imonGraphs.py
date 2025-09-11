#!/usr/bin/env python3

# Imports --------------------------------------------------
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from pathlib import Path
import argparse

# import src.helpers as hp
import helpers.accChrg as hp
from load_config import config

irradDays = re.compile(r'(?P<start>\d{6})_(?P<stop>\d{6})')
theDate = re.compile(r"(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})")


# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-d','--date', type=str, required=True, help='Stop date (YYMMDD) of irradiation period')

args = parser.parse_args()


# Processing -----------------------------------------------
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


#reads and cleans log data file using regex to only contain timestamps and imon values
vmatches,imatches = hp.parse_log(logfl)

#creates list of timesamps and imon values from cleaned log file
timestamps = hp.timestamps(imatches)
imon_values = hp.values(imatches)

#calculates accumulated charge based on imon values and timestamps
accumulated_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')

#timestamps is reported in the format datetime.datetime[year, month, day, hour, minute, second]



                                                                                            

''' shows a graph of the current vs time '''
hp.current_vs_time(start_date, end_date, timestamps, imon_values)

''' shows graph of the accumulated charge vs time '''
hp.accCharge_vs_time(start_date, end_date,timestamps, accumulated_charge)

# hp.accCharge_per_day(start_date,end_date, timestamps, accumulated_charge)

''' insteaed of daily i need an overall accumulated charge graph and calculation'''