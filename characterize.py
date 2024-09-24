#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import src.DataFile as df
import src.DataRun as dr

flnm = 'data/GasGain/test240820_S7.txt'

mdf =  df.DataFile('b904_Sr_src')

mdf.parseDataFileText(f'./{flnm}')
mruns = mdf.getDataRuns()

mdf.printRuns()
