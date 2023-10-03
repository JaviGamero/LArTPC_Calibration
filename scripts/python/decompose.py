from random import randint
from scipy.optimize import curve_fit
import numpy as np 
import matplotlib.pyplot as plt 

class greedyDecomposition(): 
    """
    The main objective of this class is to provide the necessary to separate 
    the signal of the electron from the total one using a simple method
    
    The first method: 
    For each signal, random points t will be taken to calculate the A0 value 
    from the equation 
        A(t) = A0 * e^{-t/tau},
    where tau is the decay time.
    Once A0 value is found, it will be used to check what value of the
    real signal has a larger error, this will represent the point where the 
    electron signal began. Hence, we will substract the A(t) values from this 
    point and on, resulting in an aproximation of the electron signal.
    
    The second method is to automatically fit both A0 and tau using all the 
    points of the signal. To do this, the scipy package will be used. 
    After fitting bot parameters, the process is the same.
    """
    
    def __init__(self, t, signal, GT_signal = [], tau=1.45e+03, A0=0):
        self.A0 = A0 
        self.tau = tau # (ns)
        self.t = t
        
        self.signal = signal
        self.GT_signal = GT_signal
        self.fit_signal = []
        self.decomp_signal = []

    def decayEq(self, t, A0, tau): 
        return A0 * np.exp(-t/tau)
    
    def _calculateA0(self, t, y):
        self.A0 = y * np.exp(t/self.tau)
        
    def _getA0(self):
        return self.A0
        
    def _getDecayTime(self): 
        return self.tau
        
    def manualFit(self, n = 10):
        """ 
        This function fits manually A0 considering the decay time of the 
        scintilliation ligth slow component. In order to do it, it uses n random 
        points.
        Parameters: 
        - signal: light signal to decompose
        - t: time (x axis)
        - s: seed of the random points
        - n: number of points to use when calculating A0
        """
        
        dim = len(self.signal)
        
        idxs = [randint(0,dim-1) for i in range(n)]
        A0_list = []
        
        for idx in idxs: 
            self._calculateA0(self.t[idx], self.signal[idx])
            A0_list.append(self.A0)
        
        self.A0 = np.mean(A0_list)
        self.fit_signal = self.decayEq(self.t, self.A0, self.tau)
        
    def automaticFit(self, p0 = [80, 1.6e+03]): 
        """
        This method uses an authomatic fit from the package sicpy
        """
        
        params, _ = curve_fit(self.decayEq, self.t, self.signal, p0)
        self.A0 = params[0]
        self.tau = params[1]
        self.fit_signal = self.decayEq(self.t, self.A0, self.tau)
        
    def extractElectronSignal(self, window_size = 50):
        """
        Once the experimental signal has been fitted, we first look for the 
        point with the highest error.
        One this point is calculated, we add 'window_size' (ns) to the moment 
        where the electron signal starts since scales gradually.
        """
        
        error = self.signal-self.fit_signal
        max_error_idx = np.argmax(error)
        
        e_signal = np.array(error)
        max_error_idx = np.where(self.t >= self.t[max_error_idx]-window_size)[0][0]
        e_signal[:max_error_idx] = 0 # set values before e-windowsize to 0
        
        e_signal[np.where(e_signal<0)] = 0 # negative values to 0
        
        self.decomp_signal = e_signal
    
    def plotSignals(self):
        if (np.array(self.GT_signal).size > 0): 
            fig, axs = plt.subplots(1,2, figsize = (10,5))
            
            axs[0].plot(self.t, self.signal, c='black', label='Original')
            if (np.array(self.fit_signal).size > 0): 
                axs[0].plot(self.t, self.fit_signal, c='b', label='Fit')
            if (np.array(self.decomp_signal).size > 0): 
                axs[0].plot(self.t, self.decomp_signal, c='r', label='Electron decomposed')
            axs[0].set_xlabel('Time, t (ns)')
            axs[0].set_ylabel('# Photons')
            axs[0].legend(loc='best')   
            
            axs[1].plot(self.t, self.GT_signal, c='g', label='GT')
            if (np.array(self.decomp_signal).size > 0):
                axs[1].plot(self.t, self.decomp_signal, c='r', label='Decomposition', alpha=0.8)
            axs[1].set_xlabel('Time, t (ns)')
            axs[1].set_ylabel('# Photons')
            axs[1].legend(loc='best')
            
        else: 
            fig, axs = plt.subplots(1,1, figsize = (7,7))
            
            axs.plot(self.t, self.signal, c='g', label='Original')
            if (np.array(self.fit_signal).size > 0): 
                axs.plot(self.t, self.fit_signal, c='b', label='Fit')
            if (np.array(self.decomp_signal).size > 0): 
                axs.plot(self.t, self.decomp_signal, c='r', label='Electron decomposed')
            axs.set_xlabel('Time, t (ns)')
            axs.set_ylabel('# Photons')
            # axs.set_title('Fitting')
            axs.legend(loc='best')   

        
        plt.tight_layout()
        plt.show() 
        
class quality(): 
    def __init__(self, t, GT_signal, decomp_signal):
        self.t = t
        self.GT_signal = GT_signal
        self.decomp_signal = decomp_signal
        
    def _relativeError(self, actual, calc):    
        return np.abs(actual-calc) / actual
    
    def isElectronExctracted(self, max_time_gap = 50):
        idx_max_GT = np.argmax(self.GT_signal)
        idx_max_decomp = np.argmax(self.decomp_signal)
        
        if np.abs(self.t[idx_max_decomp] - self.t[idx_max_GT]) < max_time_gap: 
            return True
        
        else: return False
        