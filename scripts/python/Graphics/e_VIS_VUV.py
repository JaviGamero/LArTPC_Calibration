import os 
import sys 

# path to python main folder in this project
libraries = os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                                         os.pardir)) 
sys.path.append(libraries) 

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import uproot

from ETL import ETL_Techniques 
from time import time

PATH = os.getcwd() # main path
print(PATH)

ROOT = os.path.join(PATH, "data/sample_particles_v2/")
print(ROOT)

################################################################################
# General variables
################################################################################
# branches in the trees
branches_to_activate_PT = ["SimPhotonsVIS", "SimPhotonsVUV", "eventID", 
                           "TrackID"]

# PMTs ids
PMTs = np.loadtxt(os.path.join(PATH, "data/PMT_IDs.txt"))
IdPMTs_L = [int(i) for i in PMTs if (i%2 == 0)] # even PMT IDs, left (X<0)
IdPMTs_R = [int(i) for i in PMTs if (i%2 != 0)] # odd PMT IDs, right (X>0)

# PMTs information --> coated/uncoated
PMTs_info = pd.read_csv(os.path.join(os.getcwd(), "data/PMT_info.txt"), 
                        sep=' ', header=0)

PMTs_info = PMTs_info[(PMTs_info['PMT_Type'] == 'pmt_coated') | 
                      (PMTs_info['PMT_Type'] == 'pmt_uncoated')]

PMTs_info = PMTs_info[['Id', 'PMT_Type']]

# time of the ideal signal
t_path = os.path.join(PATH, 'data_preproc/LightSignal_t.csv')

t = pd.read_csv(t_path, sep=';', header=None) # t[0] is nonsense, remove it
t.set_index(0, inplace=True) # remove first time

t = np.array(t).reshape(-1).astype(np.int32)
t = t-5 # move to the edges of the bins to integrate

################################################################################
# MAIN LOOP
################################################################################
skip = os.path.join(ROOT, ".DS_Store/")
count=0
df_light = pd.DataFrame()
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
            
            mu, e = etl.GTsignalsExtraction_coatedUncoated(
                eventID = branches['eventID'], 
                signalsVIS = branches['SimPhotonsVIS'], 
                signalsVUV = branches['SimPhotonsVUV'], 
                trackID = branches['TrackID'],
                PMTs_info=PMTs_info,
                light='VUV', ######## CHANGE HERE THE COMPONENT OF LIGHT TO SEE
                n_particles=4 ,
                id=1
            )
            
            e = np.array(e)[:50,:]
            mu = np.array(mu)[:50,:]
            total = e+mu
            df_light = pd.concat([df_light, pd.DataFrame(total)], axis=0)
            print(df_light.shape)
            
            
        # finish of file
        print(f'Time spent with tree {count}: {time()-mi} ')
        # print(np.shape(list_dE), np.shape(list_integrate))
        count+=1
        
print(f'Total time processing data: {time()-m0} ')

################################################################################
# Generating index
################################################################################
idx = []
for tree in range(400):
    for event in range(1,51): 
        idx.append(str(tree)+'_'+str(event))
        
################################################################################
# Saving results
################################################################################
final = df_light.set_index(pd.Index(idx)).copy()

PATH = os.getcwd() # main path
final.to_csv(os.path.join(PATH, "data_preproc/???.csv"))