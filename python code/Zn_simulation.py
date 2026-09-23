from HesseDiagram import *
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import circulant

# This module exe. the Zn simulation

# make a signal with a symmetry    
def make_signal(n_order, ZmSubgroup_gen):
    rot_arr = []
    # x = np.linspace(0, 2 * np.pi, n_order)
    y = np.random.rand(n_order)  #np.sin(x)
    n = int(n_order/ZmSubgroup_gen)
    for j in range(1,n+1):
        r_signal = np.roll(y, j*ZmSubgroup_gen)
        rot_arr.append(r_signal)
    signal = (1/n)*sum(rot_arr)
    return signal

# Step 1: initial parameters, generate signal and moments

n = 12
rho = np.random.rand(n)
rho = rho/(sum(rho))
signal_symmetry = 5 # a divisor of n
signal = make_signal(n, signal_symmetry);
# plt.plot(signal, alpha=0.9)
# plt.show()

# exact moments
c_s = circulant(signal)
M1  = c_s.dot(rho)
M2  = c_s.dot(np.diag(rho)).dot(np.transpose(c_s))
F = np.fft.fft(np.identity(n))
M1_hat = np.fft.fft(M1);
M2_hat = F.dot(M2).dot(np.matrix(F).H) # np.abs(np.fft.fft(M2))
power_spectrum = np.diagonal(M2_hat)
# plt.imshow(np.abs(M2_hat))

# Step 2: Create Hesse diagram 
generate_hasse_diagram(n)

# Step 3: Traversing the tree to get initial weights

# TBA

# Step 4: Update the weights

# TBA

# Step 5: Conclude the results


