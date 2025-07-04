#!/usr/bin/env python3

from tkinter import W
import matplotlib.pyplot as plt
import src.helpers as hp
import matplotlib.dates as mdates
from datetime import datetime

tag = '250401'

with open(f'./plots/accCharge/totalAccCharge_{tag}.txt','r') as fl:
    lines = fl.readlines()


# Initialize variables
continuous_timestamps = []
continuous_acc_charge = []
last_accumulated_charge = 0  # Starting offset for accumulated charge

step = 100000
for i,line in enumerate(lines):
    if i%step == 0:
        print(f'\tProceessig record: {i}')
    date,time,charge = line.strip().split()
    timestamp_txt = date + ' ' + time

    timestamp = datetime.fromisoformat(timestamp_txt)

    continuous_timestamps.append(timestamp)
    continuous_acc_charge.append(float(charge))


# Print overall total charge
overall_total_charge = continuous_acc_charge[-1]
print("Overall Total Charge:", overall_total_charge)

# Plot the continuous accumulated charge
plt.figure(figsize=(10, 6))
plt.plot(continuous_timestamps, continuous_acc_charge, linewidth=3, label='Accumulated Charge (mC/cm)')
plt.axvline(datetime.fromisoformat('2025-03-17 16:50:36'), color='black', linestyle='--',alpha=0.3, label='HV Scan')
print('Plotted ok')

# # Add text in the bottom-right corner
# plt.text(0.99, 0.02, f'Max Accumulated Charge Reached: {last_accumulated_charge:.2f}',
#          horizontalalignment='right',
#          verticalalignment='bottom',
#          transform=plt.gca().transAxes,
#          fontsize=8,
#          bbox=dict(facecolor='white', alpha=0.5))

# Format x-axis for datetime values
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%b-%y'))  # Date format
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=21))  # Tick every day

# Rotate x-axis labels for better readability
# plt.gcf().autofmt_xdate()

# Set labels, title, and grid
lbl_fsize = 18
plt.xlabel('Time [Date]',fontsize=lbl_fsize, ha='right', x=1.0)
plt.ylabel('Accumulated Charge [mC/cm]',fontsize=lbl_fsize, ha='right', y=1.0)
# plt.title('MiniCSC-4',fontsize=lbl_fsize)
# plt.title('Accumulated Charge Over Time',fontsize=20)
plt.tight_layout()
plt.ylim(bottom=0)
# plt.legend()
plt.grid(True, which='both', linestyle='-', alpha=0.3)  # Ensure both major/minor grids appear

hp.add_cms_label(plt.gcf(),ax=plt.gca(),context='904 Lab',cms_fontsize=lbl_fsize,status_fontsize=lbl_fsize,context_fontsize=lbl_fsize,gap=0.07)

# Show the plot
#plt.show()
plt.savefig(f'./plots/accCharge/pretty_totalAccCharge_{tag}.png',format='png',dpi=400)
plt.savefig(f'./plots/accCharge/pretty_totalAccCharge_{tag}.pdf',format='pdf')
plt.close()
