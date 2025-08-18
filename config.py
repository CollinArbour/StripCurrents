#!/usr/bin/env python
from pathlib import Path

# Configure information for run before running

# Relative Gas Gain: --------------------------------------------------
# YYYYMMDD
measurementDay = 20250808
measurementPath = f'${measurementDay}_dyn_recupCF4'


# Irradiation: --------------------------------------------------
Period01_src1_log_files = [
    './data/Irradiation/LogFiles/241018_241024_CAENGECO2020.log',
    './data/Irradiation/LogFiles/241105_241112_CAENGECO2020.log',
    './data/Irradiation/LogFiles/241115_241118_CAENGECO2020.log',
    './data/Irradiation/LogFiles/241121_241122_CAENGECO2020.log',
    './data/Irradiation/LogFiles/241122_241129_CAENGECO2020.log',
    './data/Irradiation/LogFiles/241202_241208_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250121_250206_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250206_250214_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250221_250226_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250226_250228_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250228_250303_CAENGECO2020.log',
]

Period02_src3_log_files = [
    './data/Irradiation/LogFiles/250317_250321_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250321_250322_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250323_250327_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250328_250401_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250405_250415_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250417_250429_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250429_250505_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250505_250516_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250517_250524_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250529_250602_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250604_250606_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250606_250614_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250623_250701_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250701_250711_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250714_250729_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250730_250808_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250810_250814_CAENGECO2020.log'
]

# YYMMDD
# startDay = 250730
# stopDay = 250808
# logfl = Path(f'./data/LogFiles/{startDay}_{stopDay}_CAENGECO2020.log')

# logflnm = '250604_250606_CAENGECO2020.log'
logflnm = '250810_250814_CAENGECO2020.log'
dates = logflnm.split('/')[-1].split('_')[0:2]

logfl = Path(f'./data/LogFiles/{logflnm}')
malterfl = Path(f'./data/LogFiles/MalterCurrents/{dates[1]}_Malter.log')