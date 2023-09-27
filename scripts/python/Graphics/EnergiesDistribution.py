"""
SignalComparation.py
Author: Javier Gamero Muñoz

This file will plot the distribution of two energies together: 
- "E": energy of the primary particle 
- "dE": actual energy deposited
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
from sklearn.preprocessing import MinMaxScaler 
from collections import Counter
from ETL import ETL_Techniques as etl_tech
import matplotlib.pyplot as plt 
import seaborn as sns
from time import time

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

################################################################################
# Variables
################################################################################
t0 = time()
n_tree = 0
branches_to_activate = ["stepX", "dE", "E", "PDGcode", "eventID", "energydep"]

# PMT IDs
PMTs = np.loadtxt(os.path.join(os.getcwd(), "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

# arrays that will contain the distribution of energies
dE = []
E = []
total_dE = []

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
        
        with uproot.open(file) as rootfile: 
            tree = rootfile["opanatree/OpAnaTree"] # select the tree 
            
            # reading only useful variables for this purpose
            branches = tree.arrays(branches_to_activate, library="np") 
            for entry in range(len(branches["eventID"])): # entries loop
                
                ################################################################
                # OBTAIN DISTRIBUTION OF ENERGIES
                ################################################################
                # in the decay of the muon, there are more particles apart of 
                # the Michel electron, dismiss the rest
                etl = etl_tech(IdPMTs_L, IdPMTs_R)
                dE_e, E_e, startX_e = etl.getElectron2EnergiesAndX0(branches["stepX"][entry], 
                                                               branches["PDGcode"][entry],
                                                               branches["dE"][entry],
                                                               branches["E"][entry])
                
                # and take all the energydep by all the particles
                # totalEnergyDep = etl.getTotalEnergyDep(branches["energydep"][entry])
                
                if dE_e < 0.02: continue

                dE.append(dE_e) # MeV
                # total_dE.append(totalEnergyDep) # MeV
                E.append(E_e*1e+03) # GeV to MeV
        
        n_tree += 1
        print('Tree: ', n_tree)
        
print('Time reading all the trees: ', time()-t0)

################################################################################
# PLOT
################################################################################   
# To compare both energies, we normalise them.

E_s = np.array(E).reshape(-1,1)
dE_s = np.array(dE).reshape(-1,1)

# scaler = MinMaxScaler()
# E_s = scaler.fit_transform(E_s)
# dE_s = scaler.fit_transform(dE_s)
E_s = E_s.reshape(-1)
dE_s = dE_s.reshape(-1)

data = {'E': E_s, 'dE': dE_s}
# data = {'dE': dE_s, 'total_dE': total_dE}

fig = plt.figure(figsize=(5,5))

# common_norm such that the total area of the histogram is 1 --> normalised
sns.histplot(data=data, stat='density', common_norm=True, multiple='stack', bins=500)
# plt.xlim(0.02,0.09) # when normalised
plt.xlim(2,55)
plt.xlabel("Energy, E (MeV)")
plt.ylabel("Density") 

plt.tight_layout()
plt.show()