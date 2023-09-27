"""
idealSimple_SupDecomposition.py
Author: Javier Gamero Muñoz


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
from scipy.optimize import curve_fit
from random import random, seed

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
mu_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_decomp_mu.csv')
t_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_t.csv')

mu_signals = pd.read_csv(mu_path, sep=';', header=None)
mu_signals.set_index(0, inplace=True) # set the first column as the index of each signal

t0 = 120 # (ns)
t = pd.read_csv(t_path, sep=';', header=None)
t_idx = np.where(t>t0)[1]   
t = np.array(t.iloc[0, t_idx]).reshape(-1) 

################################################################################
# MAIN
################################################################################

seed(2023)
# x = np.arange(0,1,0.01)
# y = [expDecay(i, 2, 4) + random() for i in x]

# params, _ = curve_fit(expDecay, x, y)
# print(params)
# y_fit = [expDecay(i, params[0], params[1]) for i in x]

# plt.plot(x, y, label = 'original')
# plt.plot(x, y_fit, label = 'fit')
# plt.legend(loc='best')
# plt.show()

for idx in mu_signals.index: 
    signal = mu_signals.loc[idx, t_idx].copy()
    signal = np.array(signal).reshape(-1)
    
    params, _ = curve_fit(expDecay, t, signal, [1.6e+03, 80])
    print(params)
    signal_fit = [expDecay(i, params[0], params[1]) for i in t]
    
    plt.plot(t, signal, label = 'original')
    plt.plot(t, signal_fit, label = 'fit')
    plt.ylim(0, np.max(signal))
    plt.legend(loc='best')
    
    plt.show()