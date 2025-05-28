#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import src.helpers as hp

'''
current_vs_time() takes in a date as a string in the format YEAR-MONTH-DAY
you can also put in a time if needed 

the first FULL day of data was taken on 2024-10-19
'''

# put in the dates which you need to have the data for 
log_file = './data/LogFiles/250417_250429_CAENGECO2020.log'
start_date = '2025-04-17'
end_date = '2025-04-29'


#reads and cleans log data file using regex to only contain timestamps and imon values
matches = hp.parse_log(log_file)

#creates list of timesamps and imon values from cleaned log file
timestamps = hp.timestamps(matches)
imon_values = hp.imon_values(matches)

#calculates accumulated charge based on imon values and timestamps
accumulated_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')

#timestamps is reported in the format datetime.datetime[year, month, day, hour, minute, second]



                                                                                            

''' shows a graph of the current vs time '''
hp.current_vs_time(start_date, end_date, timestamps, imon_values)

''' shows graph of the accumulated charge vs time '''
hp.accCharge_vs_time(start_date, end_date,timestamps, accumulated_charge)

# hp.accCharge_per_day(start_date,end_date, timestamps, accumulated_charge)

''' insteaed of daily i need an overall accumulated charge graph and calculation'''