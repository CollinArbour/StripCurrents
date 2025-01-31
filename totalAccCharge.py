import matplotlib.pyplot as plt
import src.helpers as hp
import matplotlib.dates as mdates

# Example Usage
log_files = [
    './data/LogFiles/241018_241024_CAENGECO2020.log',
    './data/LogFiles/241105_241112_CAENGECO2020.log',
    './data/LogFiles/241115_241118_CAENGECO2020.log',
    './data/LogFiles/241121_241129_CAENGECO2020.log'
]

# Initialize variables
continuous_timestamps = []
continuous_acc_charge = []
last_accumulated_charge = 0  # Starting offset for accumulated charge

for log in log_files:
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
plt.show()
