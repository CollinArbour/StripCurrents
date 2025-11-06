#!/usr/bin/env python
from pathlib import Path

# Chamber Specs
stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

# Configure information for run before running
stripCriteria = {'hv':3600}

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
    './data/Irradiation/LogFiles/250810_250814_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250815_250821_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250821_250827_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250827_250905_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250913_250919_CAENGECO2020.log',
    './data/Irradiation/LogFiles/250926_250930_CAENGECO2020.log',
    './data/Irradiation/LogFiles/251004_251013_CAENGECO2020.log'
]

measurementDays = [
    '20240827_dyn_recupCF4',
    # '20240916_dyn_recupCF4', # No source measurements
    '20241025_dyn_recupCF4',
    '20241031_dyn_recupCF4',
    '20241129_dyn_recupCF4',
    # '20241202_dyn_recupCF4', # No source measurements
    '20250121_dyn_recupCF4',
    '20250206_dyn_recupCF4',
    '20250221_dyn_recupCF4',
    '20250312_dyn_recupCF4',
    '20250321_dyn_recupCF4',
    '20250327_dyn_recupCF4',
    '20250401_dyn_recupCF4',
    '20250415_dyn_recupCF4',
    '20250429_dyn_recupCF4',
    '20250505_dyn_recupCF4',
    '20250516_dyn_recupCF4',
    '20250606_dyn_recupCF4',
    '20250617_dyn_recupCF4',
    '20250701_dyn_recupCF4',
    '20250711_dyn_recupCF4',
    '20250729_dyn_recupCF4',
    '20250808_dyn_recupCF4',
    '20250814_dyn_recupCF4',
    '20250821_dyn_recupCF4',
    '20250827_dyn_recupCF4',
    '20250912_dyn_recupCF4',
    '20250923_dyn_recupCF4', 
    '20251003_dyn_recupCF4', 
    '20251013_dyn_recupCF4', 
    '20251030_dyn_recupCF4', 
]

MalterChecks = [
    ('250321_Malter.log',146.60),
    ('250606_Malter.log',329.37),
    ('250729_Malter.log',427.50),
    ('250919_Malter.log',540.85),
    ('250930_Malter.log',550.69),
    ('251013_Malter.log',575.84),
]

## Dark Rates Of Interest(OI):
ratesOI = [
    '0ALCT',
    '13CLCT',
    '32TMB'
]


# ------------------------------
# Special Processing Exceptions:
# ------------------------------

# 2025-08-27
#   Dark Rate scan was done at a different time than than the Source measurements
#   Need to use run05_hv3600_dark.txt with Dark Rate Scan
#   Need to use run13_hv3600_dark.txt with Processing fits


# ----------------------------------------------------------------------------------------------------
# TO DELETE/PRUNE:
# ----------------------------------------------------------------------------------------------------

# logflnm = '250604_250606_CAENGECO2020.log'
logflnm = '250810_250814_CAENGECO2020.log'
dates = logflnm.split('/')[-1].split('_')[0:2]

logfl = Path(f'./data/LogFiles/{logflnm}')
malterfl = Path(f'./data/LogFiles/MalterCurrents/{dates[1]}_Malter.log')