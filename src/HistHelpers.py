#!/opt/anaconda3/bin/python

import numpy as np
# import matplotlib.pyplot as plt
# import uproot

def binCenters(bins):
    return 0.5 * (bins[1:] + bins[:-1])

def rebin(hist,rebin_factor):
    if len(hist[0]) % rebin_factor != 0:
        raise ValueError('Number of bins must be divisible by the rebin factor ')
    
    new_counts = np.add.reduceat(hist[0],np.arange(0,len(hist[0]),rebin_factor))
    new_bins = hist[1][::rebin_factor]
    
    return np.histogram(new_bins[:-1],new_bins,weights=new_counts)

def mGauss(x,a,x0,sig):
    return a*np.exp(-(x-x0)**2/(2*sig**2))