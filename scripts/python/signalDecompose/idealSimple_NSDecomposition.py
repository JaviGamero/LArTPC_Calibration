"""
idealSimpleDecomposition.py
Author: Javier Gamero Muñoz

This script will use a simple method to decompose the total signal and extract 
the one from the electron.

First, we will fit the exponential decay equation,  
    A(t) = A0 * e^{-gamma t},
to the slow component time serie.

After that, with the fit done for each signal, we calculate what point has the 
larger error, it will corresponds with the point where the electron decay begins
From this point and on, we remove A(t), resulting in the electron decay.
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

from decompose import idealDecompositionGreedyNS as decompose
from random import seed

################################################################################
# Functions
################################################################################

expDecay = lambda t, tau, A0: A0 * np.exp(-t/tau) # tau in (ns)

################################################################################
# Data preprocessed and variables
################################################################################
# remember that time series are loaded in each line
data_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_VUVplusVIS.csv')
t_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_t.csv')

signals = pd.read_csv(data_path, sep=';', header=None)
signals.set_index(0, inplace=True) # set the first column as the index of each signal

t0 = 200 # (ns), moment to start considering the slow component, EXPERIMENTAL
t = pd.read_csv(t_path, sep=';', header=None)
t_idx = np.where(t>t0)[1]   
t = np.array(t.iloc[0, t_idx]).reshape(-1) 

tau_slow = 1.6e+03 # decay time of the slow component

################################################################################
# Adventure begins
################################################################################
seed(2023) # results reconstruction available

for idx in signals.index: 
    signal = signals.loc[idx, t_idx].copy()
    signal = np.array(signal).reshape(-1)
    
    
    model = decompose()
    model.manualFit(signal, t, n = 300) # fit A0
    # model.automaticFit(expDecay, signal, t)
    
    Aslow = [expDecay(i, tau_slow, model.A0) for i in t] 
    
    e_signal = model.extractElectronSignal(signal, t)
    
    plt.plot(t, signal, c='g', label='Original')
    plt.plot(t, Aslow, c='b', label='Manual fit')

    plt.plot(t, e_signal, c='r', label='Electron decomposed')
    plt.legend(loc='best')
    plt.show()
    
    # TRY FITTING TAU AND SEE THE TAU 