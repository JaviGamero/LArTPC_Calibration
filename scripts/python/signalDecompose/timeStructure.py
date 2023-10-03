""" 
timeStructure.py
Author: Javier Gamero

This script will check that the time structure is preserved during the 
deconvolution process.
"""

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

from decompose import greedyDecomposition as decompose, quality
from random import seed
from time import time

################################################################################
# Functions
################################################################################

# empty 

################################################################################
# Data preprocessed and variables
################################################################################
# DECONVOLVED
m0 = time()
PATH = os.path.join(os.getcwd(), 'data_preproc/DeconvolvedSignal.csv')
decon_signals = pd.read_csv(PATH, sep=';', header=None)
decon_signals.set_index(0,inplace=True) # set the first column as the index of each signal
print("Time spent loading deconvolved data: ", time()-m0)

PATH = os.path.join(os.getcwd(), 'data_preproc/DeconvolvedSignal_t.csv')
decon_t = pd.read_csv(PATH, sep=';', header=None)
decon_t.set_index(0,inplace=True)

# IDEAL
PATH = os.path.join(os.getcwd(), 'data_preproc/LightSignal_total.csv')
ideal_signals = pd.read_csv(PATH, sep=';', header=None)
ideal_signals.set_index(0, inplace=True) 

PATH = os.path.join(os.getcwd(), 'data_preproc/LightSignal_t.csv')
ideal_t = pd.read_csv(PATH, sep=';', header=None)
ideal_t.set_index(0,inplace=True)

# Cut and stay with slow component
t0 = 200 # (ns), moment to start considering the slow component, EXPERIMENTAL
i_t_idx = np.where(ideal_t>t0)[1]
d_t_idx = np.where(decon_t>t0)[1]
ideal_t = np.array(ideal_t.iloc[0, i_t_idx]).reshape(-1) #1D array from t0 and on
decon_t = np.array(decon_t.iloc[0, d_t_idx]).reshape(-1) #1D array from t0 and on

################################################################################
# Check time structure preserves after digitalization and deconvolution
################################################################################
seed(2023) # results reconstruction available
r = [] # list of relative errors

q = quality([], [], [])
for i in range(np.shape(ideal_signals)[0]):
    print(i)
    i_signal = np.array(ideal_signals.iloc[i, i_t_idx]).reshape(-1)
    d_signal = np.array(decon_signals.iloc[i, d_t_idx]).reshape(-1)
    
    i_model = decompose(ideal_t, i_signal) 
    i_model.automaticFit() # fit the signal 
    
    d_model = decompose(decon_t, d_signal) 
    d_model.automaticFit() # fit the signal 
    
    i_tau = i_model.tau
    d_tau = d_model.tau
    
    r.append(q._relativeError(i_tau, d_tau)*100)
    
print('Mean relative error: {0}%'.format(np.mean(r))) # 7.18%