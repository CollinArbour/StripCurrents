import numpy as np
import matplotlib.pyplot as plt
from .DataRun import *

class DataFile:
    '''
    This class creates functionality for opening and parsing data files, and interacting with the run data and meta data.
    specific info about each run is stored within dataRun objects, which are held in the self.dataRuns list
    '''
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
    
    def getFileSrc(self):
        return self.dataRuns[3].getSrc()

    def filterRuns(self,use):
        '''
        Filter runs in data file by cerain conditions

        Arguments:
            @use : dict { attribute : criteria }
                e.g. {'strip':2, 'hv':3600}

        To do:
            [ ] Handle run duplicates
        '''
        to_keep = []
        to_exclude = []

        for run in self.dataRuns:
            exclude = False
            for atr in use.keys():
                if atr == 'strip':
                    if run.getStrip() != use[atr]:
                        exclude = True
                        break
                if atr == 'hv':
                    if run.getHV() != use[atr]:
                        exclude = True
                        break
            if exclude:
                to_exclude.append(run)
            else:
                to_keep.append(run)
        
        self.dataRuns = to_keep
        self.excludedRuns = to_exclude

    def sortDataRuns(self,sort):
        '''
        Sorts Data runs in a useful manner, currently supporting strip or hv

        To Do:
            [ ] Take arguments that describe condition of sorting preferred
            [ ] Throw an actual Exception if sort not supproted, don't just print
            [ ] Create a data member sort used
            [ ] Reduce coding redundancies (move sorting out of methods?)
        '''
        #Checks for valid sorting option
        Sorts = ['strip', 'hv']
        if sort not in Sorts:
            print(f'Sort option "{sort}" not supported')
            return -1

        wSrc = []
        woSrc = []

        # Seperate data runs depending on if run uses a radioactive source
        for run in self.dataRuns:
            src = run.getSrc().lower()
            if src == 'nosrc':
                woSrc.append(run)
            elif 'sr' in src or 'cd' in src:
                wSrc.append(run)

        # Sorting by Strips
        if sort == 'strip':
            if wSrc:
                wSrc_strips = []
                for run in wSrc:
                    wSrc_strips.append(run.getStrip())
                wSrc_idxs = np.argsort(wSrc_strips)         #array of indexes that the values will be taken from
                self.srcRuns = np.array(wSrc)[wSrc_idxs]    #takes the index array and sorts accordingly
            if woSrc:
                woSrc_strips = []        
                for run in woSrc:
                    woSrc_strips.append(run.getStrip())
                woSrc_idxs = np.argsort(woSrc_strips)
                self.darkRuns = np.array(woSrc)[woSrc_idxs]

        # Sorting by HV
        if sort == 'hv':
            if wSrc:
                wSrc_hvs = []
                for run in wSrc:
                    wSrc_hvs.append(run.getHV())
                wSrc_idxs = np.argsort(wSrc_hvs)
                self.srcRuns = np.array(wSrc)[wSrc_idxs]
            if woSrc:
                woSrc_hvs = []
                for run in woSrc:
                    woSrc_hvs.append(run.getHV())
                woSrc_idxs = np.argsort(woSrc_hvs)
                self.darkRuns = np.array(woSrc)[woSrc_idxs]

    def describe(self):
        print(self.name)
        if (self.srcRuns is not None):
            print('\tWith Source:')
            for run in self.srcRuns:
                print(f'\t\tRun: {run.getName()}')
                print(f'\t\t\tHV: {run.getHV()} V \t Strip: {run.getStrip()}')
                print(f'\t\t\tSrc: {run.getSrc()} \t Hole: {run.getHole()}')
        if (self.darkRuns is not None):
            print('\tWith Out Source:')
            for run in self.darkRuns:
                print(f'\t\tRun: {run.getName()}')
                print(f'\t\t\tHV: {run.getHV()} V \t Strip: {run.getStrip()}')
                print(f'\t\t\tSrc: {run.getSrc()} \t Hole: {run.getHole()}')

    def getHVScan(self,src=True):
        '''
        Returns 3 numpy arrays of hvs, avg current, and standard error

        Notes:
            Must be ran after sortDataRuns?
        '''
        hvs = []
        avgI = []
        stderr = []

        mruns = self.srcRuns if  src else self.darkRuns
        for run in mruns:
            hvs.append(run.getHV())
            avgI.append(run.getAvgCur())
            stderr.append(run.getAvgStdErr())

        return np.array(hvs),np.array(avgI),np.array(stderr)

    def getStripScan(self):
        '''
        Returns three arrays of strips, avg current, and standard error
        '''
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
        '''
            Creates dataRun objects, processess the data and meta data, and adds all dataRun objs to the self.dataRuns list
            
            Arguments:
                @mMdata: list of meta data
                @mData: run data

            Notes:
                This function is a helper function for parseDataFileText()
        '''
        self.run+=1
        datRun = DataRun(f'Run0{self.run}')
        datRun.processMetaData(mMData)
        datRun.processDataRun(mData)
        self.dataRuns.append(datRun)
        
    def parseDataFileText(self,flnm):
        '''
        This method opens a file, loops through each line, and stores the run data and meta data.

        Arguments:
            @flnm: the relative path of the file that is being read 
        '''
        #opens file and creates string array of every line           
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
                # Should execute every for loop iteration except the first
                if meta_data and run_data: 
                    self.createDataRun(meta_data,np.array(run_data).T)      #.T to group same data types together
                
                meta_data = []
                run_data = []
                meta_idxs = [i+1,i+6]       #holds line number(indexes) of all meta data lines

                #skips if statements and immediately starts next iteration
                continue

            # Store Meta Data
            # Only executes if it is within the index range determined above
            if meta_idxs and i <= meta_idxs[1] and i >= meta_idxs[0]:
                meta_data.append(line.rstrip())

            # Store numerical Data
            # Only executes if it is in the indexes after the meta_data(the actual run data) 
            if meta_idxs and i > meta_idxs[1]:
                vals = line.split()
                run_data.append(vals)
                
        if meta_data and run_data: 
            self.createDataRun(meta_data,np.array(run_data).T)
