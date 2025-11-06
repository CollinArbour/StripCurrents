import re
import sys
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_log(log_file):
    ''' reads and cleans data file returning '''
    with open(log_file, 'r') as file:
        log_data = file.read()
    
    pattern = r"\[(?P<time>.*?)\]:.*par \[(?P<label>\w+)\] val \[(?P<value>.*?)\];"
    matches = re.finditer(pattern, log_data)

    vmon = []
    imon = []

    for match in matches:
        time = match.group('time')
        value = match.group('value')

        if match.group('label') == 'VMon':
            vmon.append((time,value))
        elif match.group('label') == 'IMonH':
            imon.append((time,value))

    return vmon,imon

def timestamps(matches):
    ''' returnes a datatime array which contains timestamps in the format [datetime.datetime(year, month, day, hour, minute, second), datetime.datetime()...] '''
    # timestamps = [datetime.fromisoformat(match[0]) for match in matches]
    timestamps = [datetime.datetime.strptime(match[0],'%Y-%m-%dT%H:%M:%S') for match in matches]

    return timestamps

def values(matches):
    ''' returns the CAEN current monitor readout for plotting and calculation'''
    imon_values = [float(match[1]) for match in matches]
    
    return imon_values

def format_mmss(x,pos):
    minutes = int(x//60)
    seconds = int(x % 60)
    return f'{minutes}:{seconds:02d}'

def format_date(date):
    if len(date) == 8:
        theDate = re.compile(r"(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})")
        if mdate := re.search(theDate,date):
            return '-'.join(mdate.groups())
        else:
            sys.exit(f'ERROR with date: expected YYYYMMDD, got: {date}')

def matchAccChrg(date):
    # IMPORTANT FUNCTION TO IMPLEMENT
    return 0

def accCharge_calc(timestamps, imon_values, src='src1'):
    ''' calculates accumulated chrage based on timestamps and imon values given from log file '''

    accCharge = [0]

    if src == 'src1':
        frac = (350.05*2)/1200 / 26.22         # (I_enc * 2) / I_mon_ref for src1 / Wire Length Enclosed in 2sigma (cm)
    elif src == 'src3':
        frac = (572.54*2)/2115 / 30.48         # (I_enc * 2) / I_mon_ref for src3 / Wire Length Enclosed in 2sigma (cm)

    imon_values = np.array(imon_values)
    timestamps = np.array(timestamps)

    dts = timestamps[1:] - timestamps[:-1]
    seconds = np.array([dt.total_seconds() for dt in dts])

    daccChrg = imon_values[1:]/1000 * frac * seconds
    for dac in daccChrg:
        accCharge.append(accCharge[-1]+dac)

    return accCharge


# Rething handling of these plotting functions (move to scripts? Make it deal with axis objects instead?)
# --------------------------------------------------

def current_vs_time(start_date, end_date, timestamps, imon_values):
    ''' plots the current versus time for a specified amount of time '''
    df = pd.DataFrame({'Timestamp': timestamps, 'IMon': imon_values})
    df.set_index('Timestamp', inplace=True)

    # Filter the data based on the provided date range
    df_filtered = df[start_date:end_date]

    # Resample data to reduce noise (take mean every X minutes)
    df_resampled = df_filtered.resample('5min').mean()  # Resample per minute, adjust '1T' for different intervals (e.g., '5T' for 5 minutes)

    # Plot the resampled data
    plt.plot(df_resampled.index, df_resampled['IMon'], label="IMon (μA)", marker=".", color="blue")

    max_imonh = df_resampled['IMon'].max()  # Find the maximum IMonH value
    top_limit = max_imonh * 1.2  # Increase the y-axis limit by 20%


    # Customize the plot
    plt.title(f"IMon Over Time ({start_date} to {end_date})")
    plt.xlabel("Time (Days)")
    plt.ylabel("IMon (μA)")
    #plt.grid(True)
    plt.xticks(rotation=25)
    # plt.ylim(0, top_limit)
    plt.tight_layout()
    # plt.legend()

    # Show plot
    # plt.show()
    plt.savefig(f'./plots/accCharge/IMon/{start_date.replace('-','')}_{end_date.replace("-","")}_IMon.png', format='png', dpi=400)
    plt.close()


def accCharge_vs_time(start_date, end_date, timestamps, accumulated_charge):
    ''' plots the accumulated charge over a specified amount of time '''

    df = pd.DataFrame({'Timestamp': timestamps, 'AccumulatedCharge': accumulated_charge})
    df.set_index('Timestamp', inplace=True)

    # Filter the data based on the provided date range 
    df_filtered = df[start_date:end_date]
    
    # Plot the resampled data
    plt.plot(df_filtered.index, df_filtered['AccumulatedCharge'], label = 'Accumulated Charge (mC/cm)', marker=".", color="blue")

    max_acc_charge = df_filtered['AccumulatedCharge'].max() # find the max value of accumualted charge
    top_limit = max_acc_charge * 1.2 # increase the y-axis limit by 20%

    # Customize the plot
    plt.title(f"Accumulated Charge from ({start_date} to {end_date})")
    plt.text(0.03, 0.93, f'Accumulated {max_acc_charge:.2f} mC/cm',
        transform=plt.gca().transAxes,
        fontsize=8,
        bbox=dict(facecolor='white', alpha=0.5))
    plt.xlabel("Time (Days)")
    plt.ylabel("Accumulated Charge (mC/cm)")
    #plt.grid(True)
    plt.xticks(rotation=25)
    # plt.ylim(0, top_limit)
    plt.tight_layout()
    # plt.legend()

    # Show plot
    # plt.show()
    plt.savefig(f'./plots/accCharge/IMon/{start_date.replace("-","")}_{end_date.replace("-","")}_AccumulatedCharge.png', format='png', dpi=400)
    plt.close()