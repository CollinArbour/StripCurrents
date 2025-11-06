#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

accChrg_fl = './data_processed/measurement_accChrg.txt'
ratio24_fl = './data_processed/ratio_2_4.txt'
ratio26_fl = './data_processed/ratio_2_6.txt'

# Load Accumulated Charge Values
with open(accChrg_fl,'r') as fl:
    accChrg_data = fl.readlines()

accChrg = {}
for entry in accChrg_data[1:]:
    data = entry.split()
    accChrg[data[0]] = float(data[1])


# Load up measurement ratios
with open(ratio24_fl,'r') as fl:
    ratio24_data = fl.readlines()
ratio24_data = [line.split() for line in ratio24_data]

with open(ratio26_fl,'r') as fl:
    ratio26_data = fl.readlines()
ratio26_data = [line.split() for line in ratio26_data]


acs = []
r24 = []
r26 = []

for i in range(len(ratio24_data)):
    if ratio24_data[i][0] != ratio26_data[i][0]:
        print(f'Error, data mismatch: {ratio24_data[i][0]} != {ratio26_data[i][0]}')
    
    acs.append(accChrg[ratio24_data[i][0]])
    r24.append(float(ratio24_data[i][1]))
    r26.append(float(ratio26_data[i][1]))

plt.plot(acs,r24,marker='.',label='Ratio hole2/hole4')
plt.plot(acs,r26,marker='.',label='Ratio hole2/hole6')
plt.ylim((0.95,1.15))
plt.xlabel('Accumulated Charge (mC/cm)')
plt.legend()

plt.savefig('./data_processed/relGasGainRatios.png',format='png',dpi=300)
plt.close()