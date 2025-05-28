#!/usr/bin/env python3

import matplotlib.pyplot as plt
import src.helpers as hp
import matplotlib.dates as mdates

tag = '250505'

# Example Usage
Period01_src1_log_files = [
    './data/LogFiles/241018_241024_CAENGECO2020.log',
    './data/LogFiles/241105_241112_CAENGECO2020.log',
    './data/LogFiles/241115_241118_CAENGECO2020.log',
    './data/LogFiles/241121_241129_CAENGECO2020.log',
    './data/LogFiles/241202_241208_CAENGECO2020.log',
    './data/LogFiles/250121_250206_CAENGECO2020.log',
    './data/LogFiles/250206_250214_CAENGECO2020.log',
    './data/LogFiles/250221_250226_CAENGECO2020.log',
    './data/LogFiles/250226_250228_CAENGECO2020.log',
    './data/LogFiles/250228_250303_CAENGECO2020.log',
]

Period02_src3_log_files = [
    './data/LogFiles/250317_250321_CAENGECO2020.log',
    './data/LogFiles/250321_250322_CAENGECO2020.log',
    './data/LogFiles/250323_250327_CAENGECO2020.log',
    './data/LogFiles/250328_250401_CAENGECO2020.log',
    './data/LogFiles/250405_250415_CAENGECO2020.log',
    './data/LogFiles/250417_250429_CAENGECO2020.log',
    './data/LogFiles/250429_250505_CAENGECO2020.log',
]

# Initialize variables
continuous_timestamps = []
continuous_acc_charge = []
last_accumulated_charge = 0  # Starting offset for accumulated charge

for log in Period01_src1_log_files:
    # Parse log file
    matches = hp.parse_log(log)
    timestamps = hp.timestamps(matches)  # Extract timestamps
    imon_values = hp.imon_values(matches)  # Extract current (IMON) values

    # Calculate accumulated charge for the current file
    acc_charge = hp.accCharge_calc(timestamps, imon_values)

    # Offset accumulated charge to ensure continuity
    acc_charge = [x + last_accumulated_charge for x in acc_charge]

    # Update the last accumulated charge
    last_accumulated_charge = acc_charge[-1]

    # Append data, skipping the first timestamp and charge for subsequent files
    if continuous_timestamps:
        continuous_timestamps.extend(timestamps[1:])  # Skip the first timestamp
        continuous_acc_charge.extend(acc_charge[1:])  # Skip the first charge
    else:
        continuous_timestamps.extend(timestamps)  # Include all data for the first file
        continuous_acc_charge.extend(acc_charge)

for log in Period02_src3_log_files:
    # Parse log file
    matches = hp.parse_log(log)
    timestamps = hp.timestamps(matches)  # Extract timestamps
    imon_values = hp.imon_values(matches)  # Extract current (IMON) values

    # Calculate accumulated charge for the current file
    acc_charge = hp.accCharge_calc(timestamps, imon_values,src='src3')

    # Offset accumulated charge to ensure continuity
    acc_charge = [x + last_accumulated_charge for x in acc_charge]

    # Update the last accumulated charge
    last_accumulated_charge = acc_charge[-1]

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


# Write out total datae vs accumulated charge
tot_out = open(f'./plots/accCharge/totalAccCharge_{tag}.txt','w')
for i in range(len(continuous_timestamps)):
    tot_out.write(f'{continuous_timestamps[i]} \t {continuous_acc_charge[i]}\n')

tot_out.close()


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
plt.gcf().autofmt_xdate()

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
plt.savefig(f'./plots/accCharge/totalAccCharge_{tag}.png',format='png',dpi=400)
plt.close()
