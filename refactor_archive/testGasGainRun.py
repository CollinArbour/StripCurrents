#!/opt/anaconda3/bin/python3

# import numpy as np
# import matplotlib.pyplot as plt
# import uproot

import src.helpers as hp
import src.DAQDataRun as ddr

test = ddr.DAQDataRun('test')

test.setDataFile('./data/GasGain/fitData/output_250221_Hole2_HV3600.root')
test.setChargeSpectra()

test.fitChargeSpectra()