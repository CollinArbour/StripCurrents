#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

import src.DataFile as df
import src.DataRun as dr
import src.helpers as hp

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

for i  in range(3):
    flnm = f'./data/StripScans/Src0{i+1}.txt'

    criteria = {'hv':3600}
    
    mdf = df.DataFile(f'Sr_Src0{i+1}')
    mdf.parseDataFileText(flnm)
    mdf.filterRuns(criteria)
    mdf.sortDataRuns('strip')

    mdata = mdf.getStripScan()
    
    # Convert strip numbers into x position
    x_pos = [0.5*stripWdth]
    for j in range(1,len(mdata[0])):
        x_pos.append(x_pos[-1]+(stripWdth+stripGap))

    # Convert the strip current measurements into linear current densities
    lambda_I = np.array(mdata[1]) / stripWdth   # now has Units now nA/mm
    lambda_I_err = np.array(mdata[2]) / stripWdth   # errors also now nA/mm

    #strips = [x_pos,mdata[1],mdata[2]]
    strips = [x_pos,lambda_I,lambda_I_err]

    # Fitting the curve to sum of two Gaussians  
    # Initializing parameters with rough estimate
    p0 = [12,11,77,6,24,77]
    p1, cov = curve_fit(hp.mGaussianSum,np.array(x_pos),strips[1],p0)

    # Make plot showing 2D distribution
    hp.mkHeatMap_GaussSum(40,p1,mlabel=f'Src {i+1}',save=True)

    # Make plot displaying Strip scan shape and fit
    #hp.mkScans(strips,p1,i,save=True)
