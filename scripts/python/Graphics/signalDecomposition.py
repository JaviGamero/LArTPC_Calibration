"""
signalDecomposition.py
Author: Javier Gamero Muñoz

This script plots in two figures the signals of the mu and e added and separated
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

branches_to_activate = ["SimPhotonsVIS", "SimPhotonsVUV", "eventID", "TrackID"]

# PMT IDs
PMTs = np.loadtxt(os.path.join(os.getcwd(), "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

output = open("data_preproc/???.csv", 'w')
################################################################################
# MAIN LOOP
################################################################################
ROOT = os.path.join(os.getcwd(), "data/sample_particles_v2/")
skip = os.path.join(ROOT, ".DS_Store/")
count=0


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
                    hist_e = np.histogram(signal_e, 1000, [0,10000])[0]
                    hist_mu = np.histogram(signal_mu, 1000, [0,10000])[0]
                    t = np.arange(0,10000,10)
                    
                    
                    plt.rcParams['font.size'] = str(16)
                    
                    fig, axs = plt.subplots(1,2,figsize=(10,6))
                    plt.title('idx: '+ str(count)+'_'+str(pre_ID))                    
                    axs[0].plot(t, hist_e+hist_mu, label = 'Total', color='g')
                    axs[0].set_xlabel("Time, t (ns)")
                    axs[0].set_ylabel(r"#Photons/$\Delta$t")
                    axs[0].set_title("Total")
                    axs[0].legend(loc='best')
                    
                    axs[1].plot(t, hist_mu, alpha=0.8, color="blue", label="muon")
                    axs[1].plot(t, hist_e, color="red", label="electron")
                    axs[1].set_xlabel("Time, t (ns)")
                    axs[1].set_ylabel(r"#Photons/$\Delta$t")
                    axs[1].set_title("Signal decomposed")
                    axs[1].legend(loc='best')
                    
                    plt.tight_layout()
                    plt.show()
                    
                    # fig = plt.figure(figsize=(8,5))
                    # plt.title('idx: '+ str(count)+'_'+str(pre_ID))
                    # plt.hist(signal_e + signal_mu, 1000, [0,10000], alpha=0.5, label = 'Total', color='g')
                    # plt.hist(signal_mu, 1000, [0,10000], color="blue", label="muon")
                    # plt.hist(signal_e, 1000, [0,10000], color="red", label="electron")
                    # plt.xlabel("Time, t (ns)")
                    # plt.ylabel("# Photons")
                    # plt.legend(loc='best')          
                    
                    print('idx: ', str(count)+'_'+str(pre_ID))
                    
                    # plt.tight_layout()
                    # plt.show()
                    
                    # idx = 'idx: ' + str(count)+'_'+str(pre_ID)
                    # addToCSV(output, signal_e, len(signal_e), idx)
                    
                    signal_e = []
                    signal_mu = []
                    
                for k in range(len(branches["SimPhotonsVUV"][entry])): 
                    if (branches["TrackID"][entry] == 1): 
                        signal_mu += branches["SimPhotonsVUV"][entry][k]
                        
                    else: 
                        signal_e += branches["SimPhotonsVUV"][entry][k]
                        
                pre_ID = branches["eventID"][entry]
                
                # print("Entry: ", entry, ", eventID: ", branches["eventID"][entry], 
                #       ", pre_ID: ", pre_ID)
                
            count+=1