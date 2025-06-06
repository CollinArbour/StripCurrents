#!/opt/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import uproot
from scipy.optimize import curve_fit

import src.HistHelpers as hh
# from HistHelpers import *

class DAQDataRun:
    '''
        Class for handling data runs utilizing full DAQ system
    '''
    def __init__(self,name):
        self.name = name

        # ---- Meta Data ----
        self.date = None                # Date of run
        self.accCharge = None           # Accumulated charge

        self.runNum = None              # Run number
        self.dataFile = None            # Data file name
        self.dataTTree = None           # Data TTree
        self.darkFile = None            # Dark file name
        self.darkTTree = None           # Dark TTree

        # ---- Run Conditions ----
        self.gasComp = None             # Gas composition
        self.hv = None                  # High voltage
        self.src = None                 # Source utilized
        self.hole = None                # Hole where source is

        # ---- Data ----
        self.TMBRate = None             # TMB rate
        self.stripOccupancy = None      # Strip occupancy data
        self.wireOccupancy = None       # Wire occupancy data
        self.chargeSpectra = None       # Charge spectra data

        self.fit_db = None 

    def getName(self): return self.name
    def getDate(self): return self.date
    def getAccCharge(self): return self.accCharge
    def getRunNum(self): return self.runNum
    def getDataFile(self): return self.dataFile
    def getDarkFile(self): return self.darkFile
    def getGasComp(self): return self.gasComp
    def getHV(self): return self.hv
    def getSrc(self): return self.src
    def getHole(self): return self.hole

    # Return data histograms
    def getStripOccupancy(self): return self.stripOccupancy
    def getWireOccupancy(self): return self.wireOccupancy
    def getChargeSpectra(self): return self.chargeSpectra


    def setDate(self,date): self.date = date
    def setAccCharge(self,accCharge): self.accCharge = accCharge
    def setRunNum(self,runNum): self.runNum = runNum
    def setDataFile(self,dataFile):
        self.dataFile = dataFile
        self.dataTTree = uproot.open(dataFile)
    def setDarkFile(self,darkFile):
        self.darkFile = darkFile
        self.darkTTree = uproot.open(darkFile)

    def setChargeSpectra(self):
        if not self.dataTTree:
            print('Error, must load data file first')
            return
        self.chargeSpectra = self.dataTTree['Cathode/charge/chargeL3'].to_numpy()

    def setGasComp(self,gasComp): self.gasComp = gasComp
    def setHV(self,hv): self.hv = hv
    def setSrc(self,src): self.src = src
    def setHole(self,hole): self.hole = hole

    def saveRun(self):
        '''
        Writes DAQDataRun object to file for later use
            + Perhaps make slimmed down data to save and place in data-registry
                + This way don't need to resave most condition data
                + Store just computed quantities such as peak positions
        '''
        pass
    def loadRun(self):
        '''
        Loads DAQDataRun object that was previously stored
            + Will read in data from file and reinitialize all variables of a DAQDataRun object
            + Will open data/dark files and keep open their TTree objects
        '''
        pass

    def fitChargeSpectra(self):
        '''
        Fits charge spectra
        '''
        outpath = 'tempPlots'

        to_fit = hh.rebin(self.chargeSpectra,128)
        centers = hh.binCenters(to_fit[1])

        max_bin = np.argmax(to_fit[0])
        max_val = to_fit[0][max_bin]
        max_pos = centers[max_bin]

        fitRange=750
        fitRegion = [max_pos-fitRange,max_pos+fitRange]

        p0 = [max_val,max_pos,300]
