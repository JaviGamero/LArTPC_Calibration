from random import randint, seed
import numpy as np 

class idealDecompositionGreedyNS(): 
    """
    The main objective of this class is to provide the necessary to separate 
    the signal of the electron from the total one using a Non Supervised method
    
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
    points of the signal. To do this, the scipy package would be used. 
    """
    
    def __init__(self, decayTime=1.6e+03, A0=0):
        self.A0 = A0 
        self.decayTime = decayTime # (ns)

    def decayEq(self, t, A0=None, tau=None):
        if (A0==None and tau==None):
            return self.A0 * np.exp(-t/self.decayTime)
        else: 
            return A0 * np.exp(-t/tau)
    
    def _calculateA0(self, t, y):
        self.A0 = y * np.exp(t/self.decayTime)
        
    def _getA0(self):
        return self.A0
    
    def _calculateDecayTime(self, t1, y1, t2, y2):
        self.decayTime = (t2-t1) / np.log(y1/y2)
        
    def _getDecayTime(self): 
        return self.decayTime
        
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
        
    def extractElectronSignal(self, A0, tau, signal, t):
        """
        Once the experimental signal has been fitted, we substract it the result 
        taking into account that we cannot have negative values. 
        Parameters seem to be clear. 
        """
        
        e_signal = np.array([signal[i] - self.decayEq(t[i], A0, tau) for i in range(len(t))])
        e_signal[np.where(e_signal<0)] = 0
        
        return e_signal