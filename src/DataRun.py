import numpy as np
import matplotlib.pyplot as plt

class DataRun:
    '''
        This class is used within dataFile.py to create an object for each data run.
        The constructor takes a run name as input, and utilizes methods to declare all other member variables.
    '''
    def __init__(self,name):
        self.name = name
        #defined in processMetaData()
        self.hv = None                  #high voltage
        self.src = None                 #source
        self.strip = None               #strip
        self.layer = None               #layer
        #defined in processDataRun()
        self.vals = None                #np array of values of current
        self.timeBins = None            #np array of estimated points in time of recorded data
        self.date = None                #date of run
        self.startRun = None            #starting time of run(seconds)
        self.stopRun = None             #ending time
        #defined in characterize()
        self.avgCur = None              #average of current values
        self.stdCur = None              #standard deviation of current
        self.avgStdErr = None           #standard error 
        self.nPoints = None             #number of data points(length of vals)
        #defined in removeOutliers()
        self.nRemovedPoints = None      #number of outlier points that were removed
        self.raw = None                 #bool of whether outliers have been removed(raw = before)
        self.valsRaw = None             #storage of data before removal of outliers
        self.timeBinsRaw = None         #storage of time points before removal of outliers
        
    
    def getName(self):
        return self.name
    def getHV(self):
        return self.hv
    def getSrc(self):
        return self.src
    def getHole(self):
        return self.hole
    def getStrip(self):
        return self.strip
    
    # Handling current readings
    def getVals(self):
        return self.vals
    
    def getAvgCur(self):
        return self.avgCur
    def getStdCur(self):
        return self.stdCur
    def getAvgStdErr(self):
        return self.avgStdErr
    def getnPoints(self):
        return self.nPoints
    def getTimeBins(self):
        return self.timeBins
    
    # Handle Meta Data
    def getDate(self):
        return self.date
    def getRunStart(self):
        return self.startRun
    def getRunStop(self):
        return self.stopRun
    def getDuration(self):
        rawDur = self.stopRun-self.startRun
        return (60**2)*rawDur[0] + 60*rawDur[1] + rawDur[2]
    def getnRemoved(self):
        return self.nRemovedPoints

    def characterize(self):
        self.avgCur = np.average(self.vals)
        self.stdCur = np.std(self.vals)
        self.nPoints = len(self.vals)
        self.avgStdErr = np.abs(self.avgCur/np.sqrt(self.nPoints))
        
    def removeOutliers(self,xIQR=1.5):
        '''
            Automatically called in processDataRun(),
            Removes outliers from the data by setting an upper limit and a lower limit.
            provides functionality to save the raw data, and recalculates characterize()
            Arguments:
              @xIQR: optional parameter to set the multiplier for calculating outlier values
        '''
        #check if outliers already removed
        if not self.raw:
            print('Already removed outliers')
            return
        #first and third quartiles
        p25 = np.quantile(self.vals,0.25)
        p75 = np.quantile(self.vals,0.75)

        #middle 50%
        IQR = p75-p25

        #upper and lower limits for data
        uplim = p75 + (xIQR * IQR)
        dolim = p25 - (xIQR * IQR)

        #bool array of whether values fall within the upper and lower limits
        mmask = (self.vals < uplim) & (self.vals > dolim)

        #(num of points before) - (num of points after)
        self.nRemovedPoints = self.nPoints - sum(mmask)
        
        # Store Raw Data
        self.valsRaw = self.vals
        self.timeBinsRaw = self.timeBins

        # only saves to val if corresponding mmask element equals 'True'
        self.vals = self.vals[mmask]
        self.timeBins = self.timeBins[mmask]
        
        #redefines characterize variables and sets whether data is raw to 'false'
        self.characterize()
        self.raw = False
    
    def revertToRaw(self):
        '''
            Reverts data to before removeOutliers()
        '''
        self.vals = self.valsRaw
        self.timeBins = self.timeBinsRaw
        
        self.characterize()
        self.raw = True
    
    def processMetaData(self,mMData):
        '''
        Process the meta-data describing the Run conditions

        Arguments:
            @mMData : Array of strings corresponding to meta-data lines in text file

        Notes:
            meta-data should be properly formatted with underscores, but for now a safer more verbose processing is implemented
        '''
        #grabs meta data and splits it into a list
        #Ex:['L1', 'S1', 'HV_3600', 'Imon_1.2', 'src_Sr', 'h_02']
        meta = mMData[1] + ':' + mMData[2]
        metaspl = meta.split(':')
        
        #try-except blocks to catch improper formatting
        try:
            self.layer = int(metaspl[0].split('_')[-1])

        except:
            self.layer = int(metaspl[0][1:])
        
        try:
            self.strip = int(metaspl[1].split('_')[-1])
        except:
            self.strip = int(metaspl[1][1:])
            
        if '_' in metaspl[2]:
            self.hv = int(metaspl[2].split('_')[-1])
        else:
            self.hv = int(metaspl[2][2:])
        
        #check for source, assign src and hole accordingly
        if 'no' not in metaspl[4].lower():
            self.src = metaspl[4].split('_')[-1]
            self.hole = metaspl[5]
        else:
            self.src = 'NoSrc'
            self.hole = 'N/A'
    
    
    def processDataRun(self,mData):
        '''
            Takes list of processed data and assigns member variables respective values,
            also using helper function characterize()

            Arguments:
                @mData: list containing data for each value

            Notes:
                mData parameter should be as such:
                    -mData[0] is list of dates when value is recorded
                    -mData[1] is list of times when value is recorded
                    -mData[2] is list of actual values of current   
        '''
        
        # Store Current readings
        self.vals = np.array(mData[2],dtype=float)
        self.raw = True
        self.characterize()
        
        # Set timing information
        self.date = mData[0][0]
        self.startRun = np.array(mData[1][0].split(':'),dtype=int)
        self.stopRun = np.array(mData[1][-1].split(':'),dtype=int)
        
        # np array of estimated points in time of recorded data
        self.timeBins = np.linspace(0,self.getDuration(),num=len(self.vals))
        
        self.removeOutliers()

    
    
    def plotDataRun(self,save=False):
        '''
            Plots the current values over time.

            Arguments:
            @save: optional parameter to save the plot as a png within this directory

        '''
        # Determine where to place label in plot
        min_current = np.min(self.vals)
        max_current = np.max(self.vals)
        dif = max_current-min_current
        yloc = 0.85*dif + min_current
        
        # Create a title
        title = f'Current vs. time, Strip {self.strip}, HV {self.hv}V, Source {self.src} in hole {self.hole}'
        
        plt.plot(self.timeBins,self.vals)
        plt.axhline(self.avgCur,linestyle='--',color='black')
        plt.ylabel('Current (nA)')
        plt.xlabel('Time (s)')
        plt.text(10,yloc,f'Average Current: {self.avgCur:.4f} nA\nStandard Deviation: {self.stdCur:.4f} nA',bbox=dict(facecolor='grey',alpha=0.75))
        plt.title(title)
        
        if not save:
            plt.show()
        else:
            plt.savefig(f'./HV_{self.hv}_S{self.strip}.png',dpi=400,format='png')
        plt.close()
