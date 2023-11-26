from random import randint, seed
from scipy.optimize import curve_fit
import numpy as np 
import matplotlib.pyplot as plt 
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

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
        self.max_error = 10

    def decayEq(self, t, A0, tau): 
        return A0/tau * np.exp(-t/tau) 
    
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
        
    def _addNoiseToFitSignal(self): 
        np.random.seed(2023)
        withNoise = [i+randint(-self.max_error//2,self.max_error//2) for i in self.fit_signal]
        withNoise = np.array(withNoise)
        withNoise[withNoise<0]=0
        self.fit_signal = withNoise
    
    def plotSignals(self):
        if (np.array(self.GT_signal).size > 0): 
            plt.rcParams['font.size'] = str(16) 
            fig, axs = plt.subplots(1,2, figsize = (10,6))
            
            axs[0].plot(self.t, self.signal, c='g', label='Original', alpha=0.8)
            if (np.array(self.fit_signal).size > 0): 
                axs[0].plot(self.t, self.fit_signal, c='b', label='Muon Fit')
            if (np.array(self.decomp_signal).size > 0): 
                axs[0].plot(self.t, self.decomp_signal, c='orange', label='e (decomp)', alpha=0.7)
            axs[0].set_xlabel('Time, t (ns)')
            axs[0].set_ylabel('# Photons')
            axs[0].legend(loc='best')   
            
            axs[1].plot(self.t, self.GT_signal, c='r', label='e (GT)')
            if (np.array(self.decomp_signal).size > 0):
                axs[1].plot(self.t, self.decomp_signal, c='orange', label='e (decomp)', alpha=0.7)
            axs[1].set_xlabel('Time, t (ns)')
            axs[1].set_ylabel('# Photons')
            axs[1].legend(loc='best')
            
        else: 
            plt.rcParams['font.size'] = str(16) 
            fig, axs = plt.subplots(1,1, figsize = (8,6))
            
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
    def __init__(self, t, GT_signal=[], decomp_signal=[]):
        self.t = t
        self.GT_signal = GT_signal
        self.decomp_signal = decomp_signal
        
    def _relativeError(self, actual, calc):    
        return np.abs(actual-calc) / actual
    
    def isElectronExctracted(self, max_time_gap = 50):
        idx_max_GT = np.argmax(self.GT_signal)
        idx_max_decomp = np.argmax(self.decomp_signal)
        
        if np.abs(self.t[idx_max_decomp] - self.t[idx_max_GT]) <= max_time_gap: 
            return True
        
        else: return False
        
    def _score_efound(self, GT, pred, max_time_gap = 50):
        idx_max_GT = np.argmax(GT)
        idx_max_decomp = np.argmax(pred)
        
        if np.abs(self.t[idx_max_decomp] - self.t[idx_max_GT]) <= max_time_gap: 
            return True
        
        else: return False
        
    def getScore(self, estimator, X, y, test_size=0.2, random_state=2023, 
                 train_result=False):
        X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                            test_size=test_size, 
                                                            random_state=random_state)
        
        estimator.fit(X_train, y_train)
        y_pred_train = estimator.predict(X_train)
        y_pred = estimator.predict(X_test)
        
        e_found=0
        for GT, pred in zip(y_test, y_pred): 
            if self._score_efound(GT, pred): e_found+=1
        r_test = e_found/y_pred.shape[0]
        
        e_found_train=0
        for GT, pred in zip(y_train, y_pred_train): 
            if self._score_efound(GT, pred): e_found_train+=1
        r_train = e_found_train/y_train.shape[0]
        
        return r_train, r_test
        # if train_result: 
        
            
        
    def cross_validate(self, estimator, X, y, k=5, test_size=0.2, 
                       verbose=1, random_state = 2023): 
        
        seed(random_state)
        r_train = [] # score in train
        r_test = [] # score in test
        
        for i in range(k): 
            tr, tst = self.getScore(estimator, X, y, test_size=test_size,
                                    random_state=random_state)
            
            r_train.append(tr)
            r_test.append(tst)
            
            if verbose==1: print('Iterations of cv: {0}/{1}'.format(i+1,k))
            
        return r_train, r_test
    
    def mse(self): 
        return mean_squared_error(self.GT_signal, self.decomp_signal)
    
class calibration(): 
    """
    The intention of this class is to give it a deconvolved signal
    and convert it into the 'ideal' domain. Pass from ADC values 
    to #photons/dt.
    The method followed to calibrate the deconvolved signal is by a 
    calibration factor calculated as the division of the integral of
    the ideal series and deconvolved series in time.
    
    Parameters: 
        - signal_dec: signal(s) deconvolved to calibrate.
        
        - signal_id: signal(s) ideal, just in case to calculate calibration factor
        
        - t: time array, in case to calculate integrals
        
        - multiple: if signals_dec has more than one time serie
        
        - C: calibration factor, if None it will be calculated
        
        - pre_C: boolean, if true precalculated C will be used
    """
    def __init__(self, signal_dec, signal_id = [], t=[], multiple=True, C=None, 
                 pre_C=False):
        self.signal_dec = np.array(signal_dec)
        self.signal_id = np.array(signal_id)
        self.t = np.array(t).reshape(-1)
        self.multiple = multiple
        if pre_C: self.C = 0.016
        else: self.C = C
            
        
    def _integrateSignal(self):
        """_summary_
        Simply calculates the integral of the time series
        """
        if self.multiple:
            r_id = []
            r_dec = []
            for serie_dec, serie_id in zip(self.signal_dec, self.signal_id): 
                r_id.append(np.trapz(serie_id, self.t))
                r_dec.append(np.trapz(serie_dec, self.t))
                
            return np.array(r_id), np.array(r_dec)

        else: 
            r_id = np.trapz(self.signal_id, self.t)
            r_dec = np.trapz(self.signal_dec, self.t)
            return r_id, r_dec
        
    def _calculate_CalibrationFactor(self, bins=300, range=[0,0.2], 
                                     return_dist=False, return_C = False):
        """
        This function calulcates the calibration factor. It needs both ideal 
        signals and convolutioned signals.
        """
        
        if not self.C:  # check C is not introduced
            
            # check ideal signal and time is introduced
            if (len(self.signal_id) > 0) or (len(self.t) > 0): 
                r_id, r_dec = self._integrateSignal() # calculate integral
                self.C_distribution = r_dec/r_id # calibration factor
                
                # to calculate peak of the Gaussian
                hist = np.histogram(self.C_distribution, bins=bins, range=range) 
                max_bin = np.argmax(hist[0])
                self.C = hist[1][max_bin]
                
                # return C distribution, C or one of them
                if return_dist and return_C: 
                    print('Calibration factor and its distribution returned')
                    return self.C, self.C_distribution
                
                elif return_C: 
                    print('Calibration factor returned')
                    return self.C
                
                elif return_dist: 
                    print('Calibration factor distribution returned')
                    return self.C_distribution
                
            else: 
                print('No ideal signal introduced to calculate Calibration Factor.')
                
        else: 
            print('Calibration Factor introduced: ', self.C)

    def calibrate(self): 
        self._calculate_CalibrationFactor()
        return self.signal_dec / self.C