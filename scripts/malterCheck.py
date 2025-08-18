#!/usr/bin/env python3

import helpers.accChrg as hp
from load_config import config


# Check file exists
assert config.malterfl.exists(), f'Missing path: {config.malterfl}'

#reads and cleans log data file using regex to only contain timestamps and imon values
vmatches,imatches = hp.parse_log(config.malterfl)

hp.plotMalter(config.dates[1],vmatches,imatches)

# #creates list of timesamps and imon values from cleaned log file
# timestamps = hp.timestamps(imatches)
# imon_values = hp.imon_values(imatches)

# #calculates accumulated charge based on imon values and timestamps
# accumulated_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')

# #timestamps is reported in the format datetime.datetime[year, month, day, hour, minute, second]



                                                                                            

# ''' shows a graph of the current vs time '''
# hp.current_vs_time(start_date, end_date, timestamps, imon_values)

# ''' shows graph of the accumulated charge vs time '''
# hp.accCharge_vs_time(start_date, end_date,timestamps, accumulated_charge)

# # hp.accCharge_per_day(start_date,end_date, timestamps, accumulated_charge)

# ''' insteaed of daily i need an overall accumulated charge graph and calculation'''