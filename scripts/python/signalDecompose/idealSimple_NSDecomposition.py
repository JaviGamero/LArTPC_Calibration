"""
idealSimple_NSDecomposition.py
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

from decompose import greedyDecomposition as decompose
from random import seed

################################################################################
# Functions
################################################################################

expDecay = lambda t, A0, tau: A0 * np.exp(-t/tau) # tau in (ns)

################################################################################
# Data preprocessed and variables
################################################################################
# remember that time series are loaded in each line
data_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_total.csv')
signals = pd.read_csv(data_path, sep=';', header=None)
signals.set_index(0, inplace=True) # set the first column as the index of each signal

t_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_t.csv')
t0 = 150 # (ns), moment to start considering the slow component, EXPERIMENTAL
t = pd.read_csv(t_path, sep=';', header=None) # t[0] is nonsensen, remove it
t_idx = np.where(t>t0)[1]   
t = np.array(t.iloc[0, t_idx]).reshape(-1) #1D array from t0 and on

e_GT_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_decomp_e.csv')
e_signals = pd.read_csv(e_GT_path, sep=';', header=None)
e_signals.set_index(0, inplace=True)
e_signals = e_signals.loc[signals.index, :] # take only those with the electron

mu_GT_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_decomp_mu.csv')
mu_signals = pd.read_csv(mu_GT_path, sep=';', header=None)
mu_signals.set_index(0, inplace=True)
mu_signals = mu_signals.loc[signals.index, :] # take only those with the electron

# data does not match between GT and data without it. Rebuild it and use it:
signals = mu_signals + e_signals

################################################################################
# Adventure begins
################################################################################
seed(2023) # results reconstruction available
nphotons_min = 5 # min photons a signal needs to have in a moment to be able of 
                 # recognising the electron

for idx in signals.index: 
    print(idx)
    signal = signals.loc[idx, t_idx].copy()
    signal = np.array(signal).reshape(-1)
    
    e_signal = np.array(e_signals.loc[idx, t_idx].copy()).reshape(-1)
    
    # check electron signal is high enough to use it
    # We do this to avoid a higher error in the quality measures
    check = np.array(np.where(e_signal > 5)) 
    if not (check.size > 0): continue
    
    model = decompose(t, signal, e_signal)
    
    # # 1st method
    # model.manualFit(n = 300) # fit A0
    # signal_manualFit = model.fit_signal
    # model.extractElectronSignal()
    # e_manualSignal = model.decomp_signal
    
    # 2nd method
    model.automaticFit()
    signal_automFit = model.fit_signal
    model.extractElectronSignal()
    e_automSignal = model.decomp_signal
    
    model.plotSignals()