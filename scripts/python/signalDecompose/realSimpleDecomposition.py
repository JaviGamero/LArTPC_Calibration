"""
realSimpleDecomposition.py
Author: Javier Gamero Muñoz

This script will use a simple method to decompose the total signal and extract 
the one from the electron.

First, we will take a random moment to calculate A0 from the expression 
    A(t) = A0 * e^{-gamma t},
corresponding to the exponential decay of a particle.

After that, with an A0 for each signal, we calculate what point has the larger 
error, it will corresponds with the point where the electron decay begins.
From this point and on, we remove A(t), resulting in the electron decay.
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
import numpy as np 
import pandas as pd
from decompose import idealDecomposition as dec
import matplotlib.pyplot as plt

################################################################################
# Functions
################################################################################

# empty by the moment

################################################################################
# Data preprocessed
################################################################################
# remember that time series are loaded in each line
data_path = os.path.join(os.getcwd(), 'data_preproc/LightSignal_VUVplusVIS.csv')

signals = pd.read_csv(data_path, sep=';', header=None)
signals.set_index(0, inplace=True) # set the first column as the index of each signal

################################################################################
# Adventure begins
################################################################################
for idx in signals.index: 
    # to work, we convert each signal into a 1D numpy array
    signal = signals.loc[idx, :].copy()
    signal = np.array(signal).reshape(-1)

    # how many random points to use... mmmm
    