from random import randint, seed
from scipy.optimize import curve_fit
import numpy as np 

class idealDecompositionGreedy(): 
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
    
    def __init__(self, tau=1.6e+03, A0=0):
        self.A0 = A0 
        self.tau = tau # (ns)

    def decayEq(self, t, A0, tau): 
        return A0 * np.exp(-t/tau)
    
    def _calculateA0(self, t, y):
        self.A0 = y * np.exp(t/self.tau)
        
    def _getA0(self):
        return self.A0
        
    def _getDecayTime(self): 
        return self.tau
        
    def manualFit(self, signal, t, n = 10):
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
        
        dim = len(signal)
        
        idxs = [randint(0,dim-1) for i in range(n)]
        A0_list = []
        
        for idx in idxs: 
            self._calculateA0(t[idx], signal[idx])
            A0_list.append(self.A0)
        
        self.A0 = np.mean(A0_list)
        
    def automaticFit(self, signal, t, p0 = [80, 1.6e+03]): 
        """
        This method uses an authomatic fit from the package sicpy
        """
        
        params, _ = curve_fit(self.decayEq, t, signal, p0)
        self.A0 = params[0]
        self.tau = params[1]
        
    def extractElectronSignal(self, signal, t):
        """
        Once the experimental signal has been fitted, we first look for the 
        point with the highest error.
        From this point on, we remove the muon signal with the fitted curve. 
        Before it, we set 0.
        """
        
        signal_fit = self.decayEq(t, self.A0, self.tau)
        error = signal-signal_fit
        max_error_idx = np.argmax(error)
        
        e_signal = np.array(error)
        e_signal[:max_error_idx] = 0 # set values before e to 0
        
        e_signal[np.where(e_signal<0)] = 0 # negative values to 0
        
        return e_signal
    
    def calculateSignalFit(self, t): 
        return [self.decayEq(i, self.A0, self.tau) for i in t]