#!/usr/bin/env python3

import numpy as np

with open('./data/Irradiation/imonvalues.txt','r') as fl:
    mtxt = fl.readlines()


imons = []

for line in mtxt:
    Imon = float(line.rstrip().split()[-1])
    if not np.isnan(Imon): 
        imons.append(Imon)
        

accChrg = [0]

frac = (350.05*2)/1200

for i in imons:
    print(i)
    temp_accChrg = (i*frac/1000)/26.22*5*60
    accChrg.append(accChrg[-1]+temp_accChrg)
    
## Temporary sanity check plotting accChrg vs index
# xs = np.arange(len(accChrg))
# plt.plot(xs,accChrg)
# plt.xlabel('Seconds')
# plt.ylabel('Accumulated Charge (mC/cm)')
# plt.show()


# Writing the values to the a text file
# Superfluous for script that plots this
with open('./accChrg_vtime.txt','w') as fl:
    for i in accChrg:
        fl.write(f'{i}\n')
