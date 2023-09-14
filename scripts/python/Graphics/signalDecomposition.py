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

# empty by the moment

################################################################################
# General constants
################################################################################
baseLine = 800 # ADC value to center data and remove noise, totally expermiental
samplingTime = 2 # ns
shiftStamp = 135. # ns 

branches_to_activate = ["SimPhotonsVIS", "SimPhotonsVUV", "eventID", "TrackID"]

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
            tree = rootfile["opanatree/OpAnaPerTrackTree"] # select the tree 
            branches = tree.arrays(branches_to_activate, library="np") 
            
            signal_e = []
            signal_mu = []
            pre_ID = 1
            
            for entry in range(len(branches["eventID"])): # entries loop
                
                if (pre_ID != branches["eventID"][entry]): 
                    fig, axs = plt.subplots(1,2,figsize=(10,6))
                    axs[0].hist(signal_e + signal_mu, 1000, [0,10000], label = 'Total', color='g')
                    axs[0].set_xlabel("Time, t (ns)")
                    axs[0].set_ylabel("# Photons")
                    axs[0].set_title("Total")
                    axs[0].legend(loc='best')
                    
                    axs[1].hist(signal_mu, 1000, [0,10000], alpha=0.8, color="blue", label="muon")
                    axs[1].hist(signal_e, 1000, [0,10000], color="red", label="electron")
                    axs[1].set_xlabel("Time, t (ns)")
                    axs[1].set_ylabel("# Photons")
                    axs[1].set_title("Signal decomposed")
                    axs[1].legend(loc='best')
                    
                    fig = plt.figure(figsize=(8,5))
                    plt.hist(signal_e + signal_mu, 1000, [0,10000], alpha=0.5, label = 'Total', color='g')
                    plt.hist(signal_mu, 1000, [0,10000], color="blue", label="muon")
                    plt.hist(signal_e, 1000, [0,10000], color="red", label="electron")
                    plt.xlabel("Time, t (ns)")
                    plt.ylabel("# Photons")
                    plt.legend(loc='best')          
                    
                    plt.tight_layout()
                    plt.show()
                    
                    signal_e = []
                    signal_mu = []
                    
                for k in range(len(branches["SimPhotonsVUV"][entry])): 
                    if (branches["TrackID"][entry] == 1): 
                        signal_mu += branches["SimPhotonsVUV"][entry][k]
                        
                    else: 
                        signal_e += branches["SimPhotonsVUV"][entry][k]
                        
                pre_ID = branches["eventID"][entry]
                
                print("Entry: ", entry, ", eventID: ", branches["eventID"][entry], 
                      ", pre_ID: ", pre_ID)