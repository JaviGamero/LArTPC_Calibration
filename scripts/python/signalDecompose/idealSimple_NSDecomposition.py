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
from sklearn.metrics import mean_squared_error

from decompose import greedyDecomposition as decompose, quality
from random import seed

################################################################################
# Functions
################################################################################

expDecay = lambda t, A0, tau: A0/tau * np.exp(-t/tau) # tau in (ns)

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

r = [] # list of relative erros
e_mse = [] #list of mse
e_found = 0 # number of electrons found
e_total = 0  # number of electrons localble

for idx in signals.index: 
    print(idx)
    signal = signals.loc[idx, t_idx].copy()
    signal = np.array(signal).reshape(-1)
    
    e_signal = np.array(e_signals.loc[idx, t_idx].copy()).reshape(-1)
    mu_signal = np.array(mu_signals.loc[idx, t_idx].copy()).reshape(-1)
    
    # check electron signal is high enough to use it
    # We do this to avoid a higher error in the quality measures
    check = np.array(np.where(e_signal > 5)) 
    if not (check.size > 0): continue
    
    model = decompose(t, signal, e_signal) # if GT added, different plot
    model.automaticFit() # fit the signal 
    model.extractElectronSignal() # substract the electron signal
    e_signal_calc = model.decomp_signal
    # model.plotSignals() # plot the results
    
    model_GT = decompose(t, mu_signal)
    model_GT.automaticFit() # calculate tau 
    model.extractElectronSignal()
    
    q = quality(t, e_signal, e_signal_calc)
    r.append(q._relativeError(model_GT.tau, model.tau)*100)
    e_mse.append(q.mse())
    
    e_correct = q.isElectronExctracted()
    e_total += 1
    if e_correct: e_found+=1
    
    
print("Mean relative error: {0}%".format(np.mean(r))) # 5% of relative error --> perfect
print("MSE: {0}".format(np.mean(e_mse))) # 3.583
print('Ratio e locable found: {0}%'.format(e_found / e_total * 100)) # 72% localized --> not bad
print('Ratio e total found: {0}%'.format(e_found / signals.shape[0] * 100)) # 32%...
