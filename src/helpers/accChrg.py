import re
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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

def accCharge_calc(timestamps, imon_values, src='src1'):
    ''' calculates accumulated chrage based on timestamps and imon values given from log file '''

    accCharge = [0]

    if src == 'src1':
        frac = (350.05*2)/1200 / 26.22         # (I_enc * 2) / I_mon_ref for src1 / Wire Length Enclosed in 2sigma (cm)
    elif src == 'src3':
        frac = (572.54*2)/2115 / 30.48         # (I_enc * 2) / I_mon_ref for src3 / Wire Length Enclosed in 2sigma (cm)

    #needed to calulate the difference in seconds between log reports
    seconds = [
        (timestamps[i] - timestamps[i -1]).total_seconds()
        for i in range(1, len(timestamps))
    ]

    for i in range(len(seconds)):
        # I_mon (microA) / 1000 = I_mon (mA) which you use to report Accumulated Charge value
        temp_accCharge = (imon_values[i] * frac / 1000) * seconds[i] 
        accCharge.append(accCharge[-1] + temp_accCharge)

    return accCharge

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
    plt.ylim(0, top_limit)
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
    plt.ylim(0, top_limit)
    plt.tight_layout()
    # plt.legend()

    # Show plot
    # plt.show()
    plt.savefig(f'./plots/accCharge/IMon/{start_date.replace("-","")}_{end_date.replace("-","")}_AccumulatedCharge.png', format='png', dpi=400)
    plt.close()


def plotMalter(date,Vmon,Imon):
    # print('hello world')
    vtime = np.array(timestamps(Vmon))
    vval = values(Vmon)
    
    itime = np.array(timestamps(Imon))
    ival = np.array(values(Imon))


    if vtime[0] < itime[0]:
        offset = vtime[0]
    else:
        offset = itime[0]
    
    vtime -= offset
    itime -= offset

    vtime = [dt.total_seconds() for dt in vtime]
    itime = [dt.total_seconds() for dt in itime]

    itime=np.array(itime)
    mask = (itime > 200) * (itime < 300)
    print(itime[mask])
    print(ival[mask])

    fig, ax = plt.subplots()

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Vmon (V)',color='tab:blue')
    ax.plot(vtime,vval,color='tab:blue',label='VMon')
    ax.tick_params(axis='x',rotation=25)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_mmss))

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