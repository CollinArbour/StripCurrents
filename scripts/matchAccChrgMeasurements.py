#!/usr/bin/env python3

import re
import glob

accChrg_fl = './data/Irradiation/accumulatedCharge.txt'

measure_dir = './data/GasGain/data_processed'

# Get accumulated charge registry
with open(accChrg_fl,'r') as fl:
    accChrg = fl.readlines()

accChrg = [entry.split() for entry in accChrg[1:]]
print(f'Lenght of accChrg: {len(accChrg)}')
print('--------------------------------------------------')

# Determined processed data to match to accumulated charge values
# Collect all directories
dirs = sorted(glob.glob(f'{measure_dir}/*recupCF4'))

# Expression for collecting date info from directory name
theDate = re.compile(r"(\d{4})(\d{2})(\d{2})_dyn_recupCF4")

# Process all directories for their dates
dates = ['-'.join(date.groups())+'\n' for mdir in dirs if (date := theDate.search(mdir))]

# with open(measures_fl,'r') as fl:
#     measures = fl.readlines()

# --------------------------------------------------

periods = len(accChrg)
i=0 # irradation period index

outfl = open(f'{measure_dir}/measurement_accChrg.txt','w')
outfl.write(f'{"Measurement":<15}{"Total":>12}\n')

for measure in dates:
    measure = measure.strip()
    
    emplaced=False

    while not emplaced:
        # Check if measurement precedes irradiation:
        if i == 0 and measure < accChrg[i][0]:
            print(f'Measurement: {measure} occurs before begining of irradiation')
            outfl.write(f'{measure:<15}{0:>12.2f}\n')
            break

        # Would like to stop list index out of range errors here:
        if i > (periods-1):
            break
        elif i == (periods-1):
            if measure >= accChrg[i][2]:
                print(f'Found Measurement: {measure}')
                print(f'\tAfter: {accChrg[i][2]}')
                outfl.write(f'{measure:<15}{accChrg[i][4]:>12}\n')
            else:
                print(f'ERROR: {measure} occured during irradiation period: \n\t{accChrg[i]}')
            break

        # Increment only if needed
        if measure > accChrg[i+1][0]:
            i+=1
        
        # Check for erroneously kept log:
        if measure < accChrg[i][2]:
            print(f'ERROR: {measure} occured during irradiation period: \n\t{accChrg[i]}')
            outfl.write(f'Error with measurement: {measure}\n')
            break

        # Locate when measurement occurred
        # print(f'measure: {measure} in {i}: {accChrg[i]}')
        if i < (periods-1) and measure >= accChrg[i][2]:
            if i < (periods-1) and measure <= accChrg[i+1][0]:
                emplaced = True
                print(f'Found Measurement: {measure}')
                print(f'\tBetween: {accChrg[i][2]} and {accChrg[i+1][0]}')
                outfl.write(f'{measure:<15}{accChrg[i][4]:>12}\n')