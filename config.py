#!/usr/bin/env python
from pathlib import Path

# Configure information for run before running

# Relative Gas Gain: --------------------------------------------------
# YYYYMMDD
measurementDay = 20250808
measurementPath = f'${measurementDay}_dyn_recupCF4'


# Accumulated charge: --------------------------------------------------
# YYMMDD
# startDay = 250730
# stopDay = 250808
# logfl = Path(f'./data/LogFiles/{startDay}_{stopDay}_CAENGECO2020.log')

# logflnm = '250604_250606_CAENGECO2020.log'
logflnm = '250810_250814_CAENGECO2020.log'
dates = logflnm.split('/')[-1].split('_')[0:2]

logfl = Path(f'./data/LogFiles/{logflnm}')
malterfl = Path(f'./data/LogFiles/MalterCurrents/{dates[1]}_Malter.log')