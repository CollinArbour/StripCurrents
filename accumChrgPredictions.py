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

fitParams = []

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

    strips = [x_pos,lambda_I,lambda_I_err]

    # Fitting the curve to sum of two Gaussians  
    # Initializing parameters with rough estimate
    p0 = [12,11,77,6,24,77]
    p1, cov = curve_fit(hp.mGaussianSum,np.array(x_pos),lambda_I,p0)
    
    # Convert 1D fit into 2D distribution
    sig0 = p1[1]
    sig1 = p1[4]
    a0 = p1[0] / (np.sqrt(2*np.pi) * sig0) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = p1[3] / (np.sqrt(2*np.pi) * sig1) # Amplitude of second 2D Gaussian (nA/mm^2)
    mps = [a0,sig0,0,a1,sig1,0]

    fitParams.append(mps)

    r_names = ['Hole','FWHM','2sig','3sig','5sig']
    r_vals = [07.5,23.8/2,2*sig0,3*sig0,5*sig0]

    print(f'Source: \# {i+1}')
    for i,r in enumerate(r_vals):
        totChrg = hp.intRadius(mps,r)
        lwires = hp.wireLength(r)
        accChrg = hp.accumCharge(totChrg,r)

        actual_totChrg = hp.intRadiusCylindrical(mps,r)
        actual_accChrg = hp.accumCharge(actual_totChrg,r)

        print(f'\tBeam spot {r_names[i]} diameter: \t{2*r/10:.2f} cm')
        print(f'\t\tI enclosed  : \t\t{totChrg:.2f} nA')
        print(f'\t\tI actual  : \t\t{actual_totChrg:.2f} nA')
        print(f'\t\tNumber of wires : \t{len(hp.wirePlacement(r))}')
        print(f'\t\tLenght of wires : \t{lwires/10:.2f} cm')
        print(f'\t\tAccumulated Charge : \t{accChrg:.2f} (mC/cm) / day')
        print(f'\t\tAct. Accum. Charge : \t{actual_accChrg:.2f} (mC/cm) / day')
        print()

    

    # Make plot showing 2D distribution
    #hp.mkHeatMap_GaussSum(40,p1,mlabel=f'Src {i+1}',save=True)

    # Make plot displaying Strip scan shape and fit
    #hp.mkScans(strips,p1,i,save=False)
    #hp.mkScans(strips,mps,i,save=False,markers=False)
    #

# Adding calculations for arbitrary beam spot sizes
print('---------------------------------')
print('New Source considerations:')
print(f'Source 1 sig = {fitParams[0][1]/10:.2f} cm')
print(f'Source 3 sig = {fitParams[2][1]/10:.2f} cm')
print()

print('Src1 at d=2sig_src3: ')
r = 2*fitParams[2][1]
mps = fitParams[0]

lwires = hp.wireLength(r)
actual_totChrg = hp.intRadiusCylindrical(mps,r)
actual_accChrg = hp.accumCharge(actual_totChrg,r)

print(f'\tI Enclosed (actual): {actual_totChrg:.2f} nA')
print(f'\tWire Length: {lwires/10:.2f} cm')
print(f'\tAccumulated Charge: {actual_accChrg:.2f} (mC/cm) / day')
print('\n\n')

print('Src3 at d=2sig_src1: ')
r = 2*fitParams[0][1]
mps = fitParams[2]

lwires = hp.wireLength(r)
actual_totChrg = hp.intRadiusCylindrical(mps,r)
actual_accChrg = hp.accumCharge(actual_totChrg,r)

print(f'\tI Enclosed (actual): {actual_totChrg:.2f} nA')
print(f'\tWire Length: {lwires/10:.2f} cm')
print(f'\tAccumulated Charge: {actual_accChrg:.2f} (mC/cm) / day')
print('\n\n')