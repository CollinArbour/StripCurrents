#!/usr/bin/env python3
import os
import sys
from load_config import config

for measure in config.measurementDays:
    day = measure.split('_')[0]
    os.system(f'scripts/plotRates.py -d {day}')