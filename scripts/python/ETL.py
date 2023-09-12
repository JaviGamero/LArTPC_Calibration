import numpy as np 
import pandas as pd
import uproot

class ETL_Techniques:
    """
    ETL techniques to process data extracted from ROOT
    Author: Javier Gamero Muñoz
    
    Important, this techniques will be developed for this project, hence some 
    functions may not be useful in other jobs.
    
    Parameters
    ----------
    * 
    
    """
    
    def __init__(self, IdPMTs_L, IdPMTs_R, nEntries = 50, particlePDG = -11,
                 samplingTime=2, shiftStamp = 135., baseLine = 8000):
        
        self.IdPMTs_L = IdPMTs_L
        self.IdPMTs_R = IdPMTs_R
        self.nEntries = nEntries
        self.particlePDG = particlePDG
        
        
        self.samplingTime = samplingTime # depends on the detector
        self.shiftStamp = shiftStamp # experimental
        self.baseLine = baseLine # experimental 
        
    def getElectronEnergyAndX0(self, stepX, PDGcode, dE):
        for j in range(len(stepX)):
            if PDGcode[j]!=self.particlePDG: continue # take the electron
            dE_e = dE[j] 
            startX_e = stepX[j][0]
            
            return dE_e, startX_e
        
    def _calculateIdPMTs(self, X0):
        if X0 < 0: 
            self.sel_PMTsID = self.IdPMTs_L
        else: 
            self.sel_PMTsID = self.IdPMTs_R
            
    def getIdPMTs(self, startX_e):
        self._calculateIdPMTs(startX_e)
        return self.sel_PMTsID
    
    def getLightSignal(self, NPhotons):
        LightSignal = [] 
        for j in range(len(NPhotons)): 
            if j in set(self.sel_PMTsID): LightSignal+=NPhotons[j] 
        
        return LightSignal
    
    def getRawSignal(self, signalDigi, stampTime, pmt):
        # this is for one PMT
        x, y = [], []
        
        for k in range(len(signalDigi[pmt])): 
            rawADC = signalDigi[pmt][k] - self.baseLine
            tDigi = self.samplingTime*k + stampTime[pmt]*1000 - self.shiftStamp
            
            if (tDigi > -1000) and (tDigi<11000):
                x.append(tDigi)
                y.append(-rawADC)
                
        return [x, y]
                
    def getDecoSignal(self, signalDeco, stampTime, pmt):
        # this is for one PMT
        x, y = [], []
        
        for k in range(len(signalDeco[pmt])): 
            decoADC = signalDeco[pmt][k]/500
            tDeco = self.samplingTime*k + stampTime[pmt]*1000 - self.shiftStamp
            
            if (tDeco > -1000) and (tDeco<11000):
                x.append(tDeco)
                y.append(decoADC)
                
        return [x, y]