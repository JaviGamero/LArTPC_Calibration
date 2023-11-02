import numpy as np 
import pandas as pd
import uproot
import matplotlib.pyplot as plt

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
            dE_e = dE[j]  # this is the energy deposited
            startX_e = stepX[j][0]
            
            return dE_e, startX_e
        
    def getElectron2EnergiesAndX0(self, stepX, PDGcode, dE, E):
        # this function is just like the one above, except that it also returns
        # the real energy of the particle 
        for j in range(len(stepX)):
            if PDGcode[j]!=self.particlePDG: continue # take the electron
            dE_e = dE[j]  # this is the energy deposited
            startX_e = stepX[j][0]
            E_e = E[j]
            
            return dE_e, E_e, startX_e
        
    def getTotalEnergyDep(self, energydep, PDGcode):
        # this function calculates the sum of the energy deposition for a 
        # particle
        s = 0
        for j in range(len(PDGcode)): 
            if PDGcode[j] != self.particlePDG: continue # take the particle 
            s += np.sum(energydep[j])
            
        return s
        
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
        
        for k in range(len(signalDigi[pmt])-1): 
            rawADC = signalDigi[pmt][k] - self.baseLine
            tDigi = self.samplingTime*k + stampTime[pmt]*1000 - self.shiftStamp
            
            if (tDigi > -1000) and (tDigi<11000):
                x.append(tDigi)
                y.append(-rawADC)
                
        return [x, y]
                
    def getDecoSignal(self, signalDeco, stampTime, pmt):
        # this is for one PMT
        x, y = [], []
        
        for k in range(len(signalDeco[pmt])-1):
            decoADC = signalDeco[pmt][k]/500
            tDeco = self.samplingTime*k + stampTime[pmt]*1000 - self.shiftStamp
            
            if (tDeco > -100) and (tDeco<10000):
                x.append(tDeco)
                y.append(decoADC)
                
        return [x, y]
    
    def GTsignalsExtraction(self, eventID, signalsVIS, signalsVUV, trackID, 
                            n_particles = 4, id=1, count=0): 
        """
        This function is for one tree (per Track)
        It gets the info of the muon and electron separately for one tree
        """
        
        list_e, list_mu = [], []
        signal_e_VIS, signal_e_VUV = [], []
        signal_mu_VIS, signal_mu_VUV = [], []
        particle=0
        
        
        for entry in range(len(eventID)): # entries loop
            
            for k in range(len(signalsVUV[entry])): 
                if (trackID[entry] == id): 
                    signal_mu_VUV += signalsVUV[entry][k]
                    signal_mu_VIS += signalsVIS[entry][k]
                    
                else: 
                    signal_e_VUV += signalsVUV[entry][k]
                    signal_e_VIS += signalsVIS[entry][k]
                    
            particle+=1
            
            if (particle == n_particles): 
                particle=0
                
                idx = 'idx: ' + str(count) + '_' + str(eventID[entry])
                # print(idx)
                
                hist_VUV = np.histogram(signal_mu_VUV, 1000, [0,10000])
                hist_VIS = np.histogram(signal_mu_VIS, 1000, [0,10000])
                list_mu.append(hist_VUV[0] + hist_VIS[0])
                
                hist_VUV = np.histogram(signal_e_VUV, 1000, [0,10000])
                hist_VIS = np.histogram(signal_e_VIS, 1000, [0,10000])
                list_e.append(hist_VUV[0] + hist_VIS[0])
                
                # t = np.arange(0,10000,10)
                # plt.figure(idx)
                # plt.plot(t, hist_VUV[0] + hist_VIS[0])
                # plt.show()
                
                signal_e_VIS, signal_e_VUV = [], []
                signal_mu_VIS, signal_mu_VUV = [], []
            
        return list_mu, list_e