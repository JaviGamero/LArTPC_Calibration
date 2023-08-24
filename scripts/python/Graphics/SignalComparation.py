"""
SignalComparation.py
Author: Javier Gamero Muñoz

This file will plot in a figure the most important signals: 
    - Number of photons histogram (in the real world it cannot be obtain)
    - Time serie of the digitalized signal (this exists in the real world)
    - Time serie of the deconvolutioned signal (this exists in the real world)
"""

################################################################################
# Necessary for relative libraries
################################################################################
import os 
import sys 

libraries = os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), os.pardir)) 
sys.path.append(libraries) 

################################################################################
# Libraries 
################################################################################
import uproot 
import numpy as np 
import matplotlib.pyplot as plt 
from ETL import ETL_Techniques as etl_tech


################################################################################
# Functions
################################################################################

# empty by the moment

################################################################################
# General constants
################################################################################
baseLine = 800 # ADC value to center data and remove noise, totally expermiental
samplingTime = 2 # ns
shiftStamp = 135. # ns 

branches_to_activate = ["SimPhotonsLiteVUV", "SimPhotonsLiteVIS", 
                        "SignalsDigi", "SignalsDeco", "OpChDigi",
                        "StampTime", "StampTimeDeco",
                        "stepX", "dE", "PDGcode", "eventID"]

# PMT IDs
PMTs = np.loadtxt(os.path.join(os.getcwd(), "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

################################################################################
# MAIN LOOP
################################################################################
ROOT = os.path.join(os.getcwd(), "data/sample_particles_v2/")
skip = os.path.join(ROOT, ".DS_Store/")

for folder in os.listdir(ROOT): 
    PATH = os.path.join(ROOT, folder)
    PATH += "/"
    
    if (PATH == skip): continue # skip .DS_Store/ mac directory
    
    for f in os.listdir(PATH):
        file = os.path.join(PATH, f) # root file
        
        # another way to read variables
        with uproot.open(file) as rootfile: 
            tree = rootfile["opanatree/OpAnaTree"] # select the tree 
            
            branches = tree.arrays(branches_to_activate, library="np") 
            for entry in range(len(branches["eventID"])): # entries loop
                
                ################################################################
                # LIGHT
                ################################################################
                # in the decay of the muon, there are more particles apart of 
                # the Michel electron, dismiss the rest
                etl = etl_tech(IdPMTs_L, IdPMTs_R)
                dE_e, startX_e = etl.getElectronEnergyAndX0(branches["stepX"][entry], 
                                                       branches["PDGcode"][entry],
                                                       branches["dE"][entry])
                if dE_e < 0.0001: continue
                
                etl._calculateIdPMTs(startX_e)                
                LightSignal = etl.getLightSignal(branches["SimPhotonsLiteVUV"][entry])
                    
                ################################################################
                # ADC SIGNALS
                ################################################################
                # Calculate raw and deco signal (DIGITALIZED and DECONVOLUTIONED)
                x_raw, y_raw = [], []
                x_deco, y_deco = [], []
                for j in range(len(branches["OpChDigi"][entry])): # above PDs
                    if j not in etl.sel_PMTsID: continue
                    
                    # DIGITALIZED
                    signalDigi = etl.getRawSignal(branches["SignalsDigi"][entry], 
                                                  branches["StampTime"][entry], 
                                                  j)                    
                    x_raw += signalDigi[0]
                    y_raw += signalDigi[1]
                            
                    # DECONVOLUTIONED
                    signalDeco = etl.getDecoSignal(branches["SignalsDeco"][entry], 
                                                  branches["StampTimeDeco"][entry], 
                                                  j)
                    x_deco += signalDeco[0]
                    y_deco += signalDeco[1]
                
                ################################################################
                # PLOT
                ################################################################
                fig, axs = plt.subplots(1,3,figsize=(18, 8))
                axs[0].hist(LightSignal, 1000, [0,10000])
                axs[0].set_xlabel("Time, t (ns)")
                axs[0].set_ylabel("# Photons")
                axs[0].set_title("Light Signal (VUV, ideal)")
                
                axs[1].plot(x_raw, y_raw)
                axs[1].set_xlabel("Time, t (ns)")
                axs[1].set_ylabel("Digitalized signal, y_raw (ADC)")
                axs[1].set_title("Digitalized signal")
                
                axs[2].plot(x_deco, y_deco)
                axs[2].set_xlabel("Time, t (ns)")
                axs[2].set_ylabel("Deconvolutioned signal, y_deco (ADC)")
                axs[2].set_title("Deconvolutioned signal")
                
                plt.tight_layout()
                plt.show()