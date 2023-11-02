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

# path to python main folder in this project
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

def signal_hot(X, Y, n=3):
    
    y_max = [np.max(y) for y in Y] # max y in each pmt
    idxs = [j for j in np.argsort(y_max)[-n:]] #  take the n hottest pmt
    y_hot = [Y[i] for i in idxs] # 3 signals alone
    y_hot = np.sum([y for y in y_hot], axis=0) # 3 signals added
    
    return X[idxs[0]], y_hot

def addToCSV(f, y, dim, idx):
    f.write(idx + ';')
    
    for i in range(dim-1): 
        f.write(str(y[i]) + ';')
        
    f.write(str(y[dim-1]) + '\n')
    
    

################################################################################
# General constants
################################################################################
baseLine = 800 # ADC value to center data and remove noise, totally expermiental
samplingTime = 2 # ns
shiftStamp = 135. # ns 

branches_to_activate = ["SimPhotonsLiteVUV", "SimPhotonsLiteVIS", 
                        "SignalsDigi", "SignalsDeco", "OpChDigi",
                        "OpChDeco", "StampTime", "StampTimeDeco",
                        "stepX", "dE", "PDGcode", "eventID"]

# PMT IDs
PMTs = np.loadtxt(os.path.join(os.getcwd(), "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

output = open("data_preproc/???.csv", 'w')
n = 0
################################################################################
# MAIN LOOP
################################################################################
ROOT = os.path.join(os.getcwd(), "data/sample_particles_v2/")
skip = os.path.join(ROOT, ".DS_Store/")

for folder in os.listdir(ROOT): 
    PATH = os.path.join(ROOT, folder)
    PATH += "/"
    
    if (PATH == skip): continue # skip .DS_Store/ mac directory
    ntree = folder[-3:]
    
    for f in os.listdir(PATH):
        file = os.path.join(PATH, f) # root file
    
        
        # another way to read variables
        with uproot.open(file) as rootfile: 
            tree = rootfile["opanatree/OpAnaTree"] # select the tree 
            
            branches = tree.arrays(branches_to_activate, library="np") 
            for entry in range(len(branches["eventID"])): # entries loop
                
                event = branches["eventID"][entry]
                
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
                LightSignalVUV = etl.getLightSignal(branches["SimPhotonsLiteVUV"][entry])
                LightSignalVIS = etl.getLightSignal(branches["SimPhotonsLiteVIS"][entry])
                LightSignal = LightSignalVUV + LightSignalVIS # ideal signal is the sum up
                    
                ################################################################
                # ADC SIGNALS
                ################################################################
                # Calculate raw and deco signal (DIGITALIZED and DECONVOLUTIONED)
                
                X_raw, Y_raw = [], []
                # DIGITALIZED
                for j in range(len(branches["OpChDigi"][entry])): # above PDs
                    if (j not in etl.sel_PMTsID): continue

                    signalDigi = etl.getRawSignal(branches["SignalsDigi"][entry], 
                                                  branches["StampTime"][entry], 
                                                  j)
                    
                    if len(signalDigi[1]) == 0: continue
                        
                    X_raw.append(signalDigi[0])
                    Y_raw.append(signalDigi[1])
                    
                x_raw_hot, y_raw_hot = signal_hot(X_raw, Y_raw)
                            
                # DECONVOLUTIONED
                X_deco, Y_deco = [], []
                for j in range(len(branches["OpChDeco"][entry])): # above PDs
                    signalDeco = etl.getDecoSignal(branches["SignalsDeco"][entry], 
                                                  branches["StampTimeDeco"][entry], 
                                                  j)
                    
                    if len(signalDeco[1]) == 0: continue
                    
                    X_deco.append(signalDeco[0])
                    Y_deco.append(signalDeco[1])                    
                    
                x_deco_hot, y_deco_hot = signal_hot(X_deco, Y_deco, n=1)
                
                # if (len(x_deco_hot) != 5050 or len(y_deco_hot) != 5050): 
                #     print(str(n) + '_' + str(event))
                #     print(len(x_deco_hot), len(y_deco_hot))
                
                # addToCSV(output, y_deco_hot, len(y_deco_hot), 
                        #  str(n) + '_' + str(event))
                
                # addToCSV(output, x_deco_hot, len(x_deco_hot), '0')
                # input("Time array loaded. Stop.")
                    
                ################################################################
                # PLOT
                ################################################################
                print(str(n) + '_' + str(event))
                fig, axs = plt.subplots(1,3,figsize=(18, 8))
                axs[0].hist(LightSignal, 1000, [0,10000], color = 'g', label='Total')
                axs[0].set_xlabel("Time, t (ns)")
                axs[0].set_ylabel("# Photons")
                axs[0].set_title("Light Signal (ideal)")
                axs[0].legend(loc='best')
                
                axs[1].plot(x_raw_hot, y_raw_hot, color = 'g', label='Total')
                axs[1].set_xlabel("Time, t (ns)")
                axs[1].set_ylabel("Digitalized signal, y_raw (ADC)")
                axs[1].set_title("Digitalized signal")
                axs[1].legend(loc='best')
                
                axs[2].plot(x_deco_hot, y_deco_hot, color = 'g', label='Total')
                axs[2].set_xlabel("Time, t (ns)")
                axs[2].set_ylabel("Deconvolutioned signal, y_deco (ADC)")
                axs[2].set_title("Deconvolutioned signal")
                axs[2].legend(loc='best')
                
                plt.tight_layout()
                plt.show()
            
            n+=1
            print('Tree: ', n)

output.close()