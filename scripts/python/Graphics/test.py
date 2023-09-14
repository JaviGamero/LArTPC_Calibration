import numpy as np 
from collections import Counter
from sklearn.preprocessing import MinMaxScaler

E = [1,2,3,4,5,3,4,5,2,3,2,1,3,4,5]
scaler = MinMaxScaler()

a = np.array(E).reshape(-1,1)
print(a)
print(len(a))
print(a.reshape(1, len(a)))