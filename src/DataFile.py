import numpy as np
import matplotlib.pyplot as plt
from .DataRun import *

class DataFile:
    def __init__(self,name):
        self.name = name
        self.run = -1
        self.dataRuns = []
    
    def printRuns(self):
        to_return = ''
        for run in self.dataRuns:
            to_return += f'{run.getName()} '
        print(to_return)
    
    def getDataRuns(self):
        return self.dataRuns
    
    def createDataRun(self,mMData,mData):
        self.run+=1
        datRun = DataRun(f'Run0{self.run}')
        datRun.processMetaData(mMData)
        datRun.processDataRun(mData)
        self.dataRuns.append(datRun)
        
    def parseDataFileText(self,flnm):
        with open(flnm,'r') as fl:
            textData = fl.readlines()
        
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
    
