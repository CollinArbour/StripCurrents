#!/usr/bin/env python3

# Imports --------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse

import helpers.accChrg as hp
from load_config import config


# Arguments ------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-t','--tag', type=str, default='TEMP', help='Tag name to give saved output')

args = parser.parse_args()

print(f' ~~ Storing with tag: {args.tag} ~~\n')


# Processing -----------------------------------------------
# Initialize variables
continuous_timestamps = []
continuous_acc_charge = []
last_accumulated_charge = 0  # Starting offset for accumulated charge

# Create record keeping for accumulated charge at each point
logfl = open('./data/Irradiation/accumulatedCharge.txt','w')
logfl.write(f'{"Dates":<25}{"Accumulated":>15}{"Total":>15}\n')

def updateLog(log,accumulated,total,verbose=False):
    start,stop = log.split('/')[-1].split('_')[:2]
    fstart = f'20{start[:2]}-{start[2:4]}-{start[4:]}'
    fstop = f'20{stop[:2]}-{stop[2:4]}-{stop[4:]}'
    range = f'{fstart} to {fstop}'

    if verbose:
        print(f'\t{range:<25}{accumulated:>15.2f}{total:>15.2f}')
    logfl.write(f'{range:<25}{accumulated:>15.2f}{total:>15.2f}\n')


print('Processing Source1 Period') 
for log in config.Period01_src1_log_files:
    # Parse log file
    matches = hp.parse_log(log)
    timestamps = hp.timestamps(matches[1])  # Extract timestamps
    imon_values = hp.values(matches[1])  # Extract current (IMON) values

    # Calculate accumulated charge for the current file
    acc_charge = hp.accCharge_calc(timestamps, imon_values)
    accumulated = acc_charge[-1]

    # Offset accumulated charge to ensure continuity
    acc_charge = [x + last_accumulated_charge for x in acc_charge]

    # Update the last accumulated charge
    last_accumulated_charge = acc_charge[-1]

    # Log contribution to accumulatedCharge.txt
    updateLog(log,accumulated,last_accumulated_charge,verbose=True)
    
    # Append data, skipping the first timestamp and charge for subsequent files
    if continuous_timestamps:
        continuous_timestamps.extend(timestamps[1:])  # Skip the first timestamp
        continuous_acc_charge.extend(acc_charge[1:])  # Skip the first charge
    else:
        continuous_timestamps.extend(timestamps)  # Include all data for the first file
        continuous_acc_charge.extend(acc_charge)

print('Processing Source3 Period') 
for log in config.Period02_src3_log_files:
    # Parse log file
    matches = hp.parse_log(log)
    timestamps = hp.timestamps(matches[1])  # Extract timestamps
    imon_values = hp.values(matches[1])  # Extract current (IMON) values

    # Calculate accumulated charge for the current file
    acc_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')
    accumulated = acc_charge[-1]

    # Offset accumulated charge to ensure continuity
    acc_charge = [x + last_accumulated_charge for x in acc_charge]

    # Update the last accumulated charge
    last_accumulated_charge = acc_charge[-1]

    # Log contribution to accumulatedCharge.txt
    updateLog(log,accumulated,last_accumulated_charge,verbose=True)

    # Append data, skipping the first timestamp and charge for subsequent files
    if continuous_timestamps:
        continuous_timestamps.extend(timestamps[1:])  # Skip the first timestamp
        continuous_acc_charge.extend(acc_charge[1:])  # Skip the first charge
    else:
        continuous_timestamps.extend(timestamps)  # Include all data for the first file
        continuous_acc_charge.extend(acc_charge)



# Print overall total charge
overall_total_charge = continuous_acc_charge[-1]
print("Overall Total Charge:", overall_total_charge)

# Plot the continuous accumulated charge
plt.figure(figsize=(10, 6))
plt.plot(continuous_timestamps, continuous_acc_charge, label='Accumulated Charge (mC/cm)')

# Add text in the bottom-right corner
plt.text(0.99, 0.02, f'Max Accumulated Charge Reached: {last_accumulated_charge:.2f}',
         horizontalalignment='right',
         verticalalignment='bottom',
         transform=plt.gca().transAxes,
         fontsize=8,
         bbox=dict(facecolor='white', alpha=0.5))

# Format x-axis for datetime values
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))  # Date format
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))  # Tick every day

# Rotate x-axis labels for better readability
plt.gcf().autofmt_xdate(rotation=50)

# Set labels, title, and grid
plt.xlabel('Time (Days)')
plt.ylabel('Accumulated Charge (mC/cm)')
plt.title('Accumulated Charge Over Time')
plt.tight_layout()
plt.ylim(bottom=0)
plt.legend()
plt.grid(True, which='both', linestyle='-', alpha=0.8)  # Ensure both major/minor grids appear

# Show the plot
#plt.show()
plt.savefig(f'./plots/accCharge/totalAccCharge_{args.tag}.png',format='png',dpi=400)
plt.close()

logfl.close()