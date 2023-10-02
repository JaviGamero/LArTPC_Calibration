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
m0 = time()
decon_path = os.path.join(os.getcwd(), 'data_preproc/DeconvolvedSignal.csv')
decon_signals = pd.read_csv(decon_path, sep=';', header=None)
print("Time spent loading data: ", time()-m0)
print(decon_signals.head())
print(decon_signals.shape)