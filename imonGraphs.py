import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import src.helpers as hp


'''
this is where we are going to call the imon graphs for right now

current_vs_time() takes in a date as a string in the format YEAR-MONTH-DAY
you can also put in a time if needed 

the first FULL day of data was taken on 2024-10-19

'''

# put in the dates which you need to have the data for 
start_date = '2024-10-18'
end_date = '2024-10-24'

#CAEN log reading and parsing
with open('./data/LogFiles/CAENGECO2020.log', 'r') as file:
        log_data = file.read()

#cleans log data file using regex to only contain timestamps and imon values
matches = hp.parse_log(log_data)

#creates list of timesamps and imon values from cleaned log file
timestamps = hp.timestamps(matches)
imon_values = hp.imon_values(matches)










hp.current_vs_time(start_date, end_date, timestamps, imon_values)


#hp.accCharge_vs_time(start_date, end_date)