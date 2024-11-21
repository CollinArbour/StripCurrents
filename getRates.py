import numpy as np
import src.helpers as hp

#defines the specific TMB dump you want to process
tmbbasedir = './data/GasGain/TMBDumps'
run_dir = '240826_HV3600'

#removes extranious files that may not be tmb dumps and returns array of files sorted into data and dark
mfd_runs,mfb_runs = hp.mruns(tmbbasedir,run_dir) 

#extracts ALCT, CLCT, and TMB data from each dump file
mfd_rates = hp.processTMBDumps(tmbbasedir,run_dir,mfd_runs)
mfb_rates = hp.processTMBDumps(tmbbasedir,run_dir,mfb_runs)

#print(f'Run: {mfd_rates[0][0]}\t TMB Rate: {mfd_rates[1][0]} Hz')
#print(f'Run: {mfb_rates[0][0]}\t TMB Rate: {mfb_rates[1][0]} Hz')
#print(f'Ratio ( bkg/(s+b) ): \t\t {mfb_rates[1][0][-1]/mfd_rates[1][0][-1]}')

#for i in range(len(mfb_rates[1])):
 #   print(f'{mfd_rates[1][i]}') 

print(f'Data Runs for {run_dir}')
for i in range(len(mfd_runs)):
    print(f'{i}\tRun: {mfd_rates[0][i]}\t TMB Rate: {mfd_rates[1][i][-1]} Hz')


print(f'Dark Runs for {run_dir}')
for i in range(len(mfb_runs)):
    print(f'{i}\tRun: {mfb_rates[0][i]}\t TMB Rate: {mfb_rates[1][i][-1]} Hz')

'''
old code from previous iteration saving for 

#directory = '241025_firstRefMeasures'
#print(hp.mruns(tmbbase,directory))

base_dir = "./data/GasGain/TMBDumps/240916_TMB_Dumps"

parsed_data = hp.parse_directory(base_dir)

for entry in parsed_data:
    print(f'{entry}\t TMB Rate: {entry['tmb_rate']}')

if you want to report all data

for entry in parsed_data:
    print(f"File: {entry['file_path']}")
    print(f"  ALCT Rate: {entry['alct_rate']}")
    print(f"  CLCT Rate: {entry['clct_rate']}")
    print(f"  TMB Rate: {entry['tmb_rate']}")
    print()


if you just want the tmb rates (shitty dont use this )
for entry in parsed_data:
    print(f'File: {entry['file_path']}')
    print(f'  TMB Rate: {entry['tmb_rate']}')
    print()

'''
