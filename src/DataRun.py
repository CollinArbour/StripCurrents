import numpy as np
import matplotlib.pyplot as plt

class DataRun:
    def __init__(self,name):
        self.name = name
    
    def getName(self):
        return self.name
    def getHV(self):
        return self.hv
    def getSrc(self):
        return self.src
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
        self.avgStdErr = self.avgCur/np.sqrt(self.nPoints)
        
    def removeOutliers(self,xIQR=1.5):
        if not self.raw:
            print('Already removed outliers')
            return
        
        p75 = np.quantile(self.vals,0.75)
        p25 = np.quantile(self.vals,0.25)

        IQR = p75-p25

        uplim = p75 + (xIQR * IQR)
        dolim = p25 - (xIQR * IQR)

        upmask = self.vals < uplim
        domask = self.vals > dolim

        mmask = upmask * domask
        
        self.nRemovedPoints = self.nPoints - sum(mmask)
        
        # Store Raw Data
        self.valsRaw = self.vals
        self.timeBinsRaw = self.timeBins

        self.vals = self.vals[mmask]
        self.timeBins = self.timeBins[mmask]
        
        self.characterize()
        self.raw = False
    
    def revertToRaw(self):
        self.vals = self.valsRaw
        self.timeBins = self.timeBinsRaw
        
        self.characterize()
        self.raw = True
    
    def processMetaData(self,mMData):
        # Start processing Meta Data
        meta = mMData[1] + ':' + mMData[2]
        metaspl = meta.split(':')
        
        try:
            self.layer = int(metaspl[0].split('_')[-1])
        except:
            self.layer = int(metaspl[0][1:])
        
        # Get strip number, convention changed and I am employing this safety strategy
        try:
            self.strip = int(metaspl[1].split('_')[-1])
        except:
            self.strip = int(metaspl[1][1:])
            
            
        self.hv = int(metaspl[2].split('_')[-1])
        self.src = metaspl[4].split('_')[-1]
        self.hole = metaspl[5].split('_')[-1]
#         self.hole = int(metaspl[5].split('_')[-1])
    
    
    def processDataRun(self,mData):
        
        # Store Current readings
        self.vals = np.array(mData[2],dtype=float)
        self.raw = True
        self.characterize()
        
        # Set timing information
        self.date = mData[0][0]
        self.startRun = np.array(mData[1][0].split(':'),dtype=int)
        self.stopRun = np.array(mData[1][-1].split(':'),dtype=int)
        
        # Determine the time bins
        self.timeBins = np.linspace(0,self.getDuration(),num=len(self.vals))
        
        self.removeOutliers()

    
    
    def plotDataRun(self,save=False):
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
