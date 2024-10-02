import numpy as np
import matplotlib.pyplot as plt
from .DataRun import *

class DataFile:
    def __init__(self,name):        #constructor
        self.name = name
        self.run = -1
        self.dataRuns = []

        self.srcRuns = None
        self.darkRuns = None
    
    def printRuns(self):
        to_return = ''
        for run in self.dataRuns:
            to_return += f'{run.getName()} '
        print(to_return)
    
    def getDataRuns(self):
        return self.dataRuns

    def sortDataRuns(self):
        wSrc = []
        woSrc = []
        for run in self.dataRuns:
            src = run.getSrc().lower()
            if src == 'no':
                woSrc.append(run)
            elif 'sr' in src or 'cd' in src:
                wSrc.append(run)

        if wSrc:
            wSrc_strips = []
            
            for run in wSrc:
                wSrc_strips.append(run.getStrip())

            wSrc_idxs = np.argsort(wSrc_strips)
            self.srcRuns = np.array(wSrc)[wSrc_idxs]

        if woSrc:
            woSrc_strips = []        

            for run in woSrc:
                woSrc_strips.append(run.getStrip())

            woSrc_idxs = np.argsort(woSrc_strips)
            self.darkRuns = np.array(woSrc)[woSrc_idxs]

    def getStripDist(self):
        if self.srcRuns is None:
            print('No runs with source recorded')
            return
        
        strips = []
        avgI = []
        stderr = [] 

        for run in self.srcRuns:
            strips.append(run.getStrip())
            avgI.append(run.getAvgCur())
            stderr.append(run.getAvgStdErr())

        return strips,avgI,stderr
    
    def createDataRun(self,mMData,mData):
        self.run+=1
        datRun = DataRun(f'Run0{self.run}')
        datRun.processMetaData(mMData)
        datRun.processDataRun(mData)
        self.dataRuns.append(datRun)
        
    def parseDataFileText(self,flnm):           #open file(flnm)
        with open(flnm,'r') as fl:              #renamed fl
            textData = fl.readlines()           #create list of strings, each line
        
        # Initialize temp storage arrays
        meta_data = []
        run_data = []
        meta_idxs = []
        
        # Loop over lines of text data file
        for i,line in enumerate(textData):
            # Look for text marking beginning of run
            if '# v' in line:

                # If a previous data run has been recorded
                # Cache it before clearing and beginning new collection
                if meta_data and run_data: 
                    self.createDataRun(meta_data,np.array(run_data).T)
                
                meta_data = []
                run_data = []
                meta_idxs = [i+1,i+6]

                continue

            # Store Meta Data
            if meta_idxs and i <= meta_idxs[1] and i >= meta_idxs[0]:
                meta_data.append(line.rstrip())

            # Store numerical Data
            if meta_idxs and i > meta_idxs[1]:
                vals = line.split()
                run_data.append(vals)
                
        if meta_data and run_data: 
            self.createDataRun(meta_data,np.array(run_data).T)
