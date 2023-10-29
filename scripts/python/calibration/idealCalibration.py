################################################################################
# Necessary for relative libraries
################################################################################
import os 
import sys 

# path to python main folder in this project
libraries = os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                                         os.pardir)) 
sys.path.append(libraries) 

################################################################################
# Libraries 
################################################################################
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import uproot

from ETL import ETL_Techniques 
from decompose import calibration
from random import seed
from time import time

################################################################################
# Functions
################################################################################


################################################################################
# General variables
################################################################################
# branches in the trees
branches_to_activate_PT = ["SimPhotonsVIS", "SimPhotonsVUV", "eventID", 
                           "TrackID"]
branches_to_activate_total = ["stepX", "dE", "PDGcode", "eventID"]

# PMTs ids
PMTs = np.loadtxt(os.path.join(os.getcwd(), "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

# time of the ideal signal
PATH = os.getcwd()
t_path = os.path.join(PATH, 'data_preproc/LightSignal_t.csv')

t = pd.read_csv(t_path, sep=';', header=None) # t[0] is nonsense, remove it
t.set_index(0, inplace=True) # remove first time

t = np.array(t).reshape(-1).astype(np.int32)
t = t-5 # move to the edges of the bins to integrate
################################################################################
# MAIN LOOP
################################################################################

ROOT = os.path.join(os.getcwd(), "data/sample_particles_v2/")
skip = os.path.join(ROOT, ".DS_Store/")
count=0
list_dE, list_integrate = [], []
m0 = time()

etl = ETL_Techniques(IdPMTs_L, IdPMTs_R)

for folder in os.listdir(ROOT): 
    PATH = os.path.join(ROOT, folder)
    PATH += "/"
    
    if (PATH == skip): continue # skip .DS_Store/ mac directory
    
    for f in os.listdir(PATH):
        file = os.path.join(PATH, f) # root file
        mi = time()
        
        # open trees file:
        with uproot.open(file) as rootfile: 
            # first tree: per particle
            tree = rootfile["opanatree/OpAnaPerTrackTree"]
            branches = tree.arrays(branches_to_activate_PT, library="np") 
            
            list_mu, list_e = etl.GTsignalsExtraction(eventID = branches['eventID'], 
                                                      signalsVIS = branches['SimPhotonsVIS'], 
                                                      signalsVUV = branches['SimPhotonsVUV'], 
                                                      trackID = branches['TrackID'], 
                                                      id=1, count=count)
            
            
            # second tree: with total signal
            tree = rootfile["opanatree/OpAnaTree"]
            branches = tree.arrays(branches_to_activate_total, library="np") 
            dE = []
            
            for entry in range(len(branches["eventID"])): # entries loop
                event = branches["eventID"][entry]
                
                ################################################################
                # LIGHT
                ################################################################
                # in the decay of the muon, there are more particles apart of 
                # the Michel electron, dismiss the rest
                etl = ETL_Techniques(IdPMTs_L, IdPMTs_R)
                dE_e, _ = etl.getElectronEnergyAndX0(branches["stepX"][entry], 
                                                       branches["PDGcode"][entry],
                                                       branches["dE"][entry])
                dE.append(dE_e)
                
            cal = calibration(list_e, t, multiple=True)
            list_integrate.append(cal._integrateSignal())
            list_dE.append(dE)
                
            
        # finish of file
        print(f'Time spent with tree {count}: {time()-mi} ')
        # print(np.shape(list_dE), np.shape(list_integrate))
        count+=1
        
print(f'Total time processing data: {time()-m0} ')
list_dE = np.array(list_dE).reshape(-1)
list_integrate = np.array(list_integrate).reshape(-1)

order = np.argsort(list_dE)
list_dE = list_dE[order]
list_integrate = list_integrate[order]

plt.figure()
plt.plot(list_dE, list_integrate, c='b')
plt.xlabel('Energy (MeV)')
plt.ylabel('Integral of the e signal')
plt.show()