import numpy as np 

class idealDecompositionNS(): 
    """
    The main objective of this class is to provide the necessary to separate 
    the signal of the electron from the total one using a Non Supervised method
    
    For each signal, a random point t will be taken to calculate the A0 value 
    from the equation 
        A(t) = A0 * e^{-gamma t},
    where gamma is the decay time
    
    Once A0 value is found, it will be used to check what value of the
    real signal has a larger error, this will represent the point where the 
    electron signal began. Hence, we will substract the A(t) values from this 
    point and on, resulting in an aproximation of the electron signal.
    """
    
    def __init__(self, decayTime, A0=0):
        self.A0 = A0
        self.decayTime = decayTime
    
    def _calculateA0(self, t, y):
        self.A0 = y * np.exp(t/self.decayTime)
        
    def _getA0(self):
        return self.A0