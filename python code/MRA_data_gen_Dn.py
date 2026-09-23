# This code generate symmetric signals and their respective 1st and 2nd moments
#
#       The case of the Dihedral group
#
# NS, September 2024

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import circulant


# adding noise according to required SNR
def add_gaussian_noise(signal, snr):
    # Calculate signal power and convert SNR from dB to linear scale
    signal_power = np.mean(signal**2)
    snr_linear = 10**(snr / 10)

    # Calculate noise power
    noise_power = signal_power / snr_linear

    # Generate Gaussian noise
    noise = np.sqrt(noise_power) * np.random.randn(len(signal))

    # Add noise to the signal
    noisy_signal = signal + noise
    return noisy_signal

# Generate synthetic data
def generate_data(num_references, signal, noise_snr, rho):
    data = []
    num_points =  np.size(signal)
    y = signal #np.sin(x)
    values = range(0,num_points)
    group_sampls = random.choices(values, rho, k=num_references);
    
    for el in range(num_references):
        shift_amount = group_sampls[el] # np.random.randint(0,num_points-1)
        y_rolled = np.roll(y, shift_amount)
        noisy_signal = add_gaussian_noise(y_rolled, noise_snr)
        data.append(noisy_signal)
    return np.array(data)

# Align data by shifting to the mean
def align_data(data):
    mean_reference = np.mean(data, axis=0)
    aligned_data = []
    for y in data:
        shift = np.argmax(np.correlate(y, mean_reference, mode='full')) - len(mean_reference) + 1
        aligned_data.append(np.roll(y, -shift))
    return np.array(aligned_data)

# Plot data
def plot_data(data, title):
    for y in data:
        plt.plot(y, alpha=0.5)
    plt.title(title)
    plt.show()
    
# make a signal with shift symmetry    
def make_signal(n, ZmSubgroup_gen):
    rot_arr = []
    y = np.random.rand(n)  #np.sin(x)
    ell = int(n/ZmSubgroup_gen)
    for j in range(1,ell+1):
        r_signal = np.roll(y, j*ZmSubgroup_gen)
        rot_arr.append(r_signal)
    signal = (1/n)*sum(rot_arr)
    return signal

def plot_discrete(vec, nameit):
    
    n = np.size(vec)
    plt.stem(range(n), vec, use_line_collection=True)
    #plt.xlabel('n')
    #plt.xticks(np.arange(LL, UL, 1))
    #plt.yticks([0, 1])
    plt.ylabel('value')
    plt.title(nameit)
    plt.show()

# OLD VERSION
# def make_signal_ref(num_points,m):
#     # y = np.arange(1, num_points + 1) # 
#     y = np.random.rand(num_points)
#     # Reflect the signal with respect to the middle
#     middle_index = len(y) // 2
#     sy = np.concatenate((y[middle_index:][::-1], y[:middle_index][::-1]))
#     #sy2 = y[:][::-1]
#     signal = (1/2)*(y+sy)
#     rot_arr = []
#     # rotating
#     n = int(num_points/m)
#     for j in range(1,n+1):
#         r_signal = np.roll(signal, j*m)
#         rot_arr.append(r_signal)
#     signal_sym = (1/n)*sum(rot_arr)
#     return signal_sym


# Generating a random signal where the symmetry subgroup is 
# D_mj = { r^n/m,r^2(n/m),..,1,r^(j+n/m)s,..,r^(j+n-n/m)s }
# m is the order, that is D_mj=~D2m
def make_sym_signal(n, m, j):
    y = np.random.randn(n)  # a random vector
    rot_arr = []                    # initializing
    # the "m" part: 
    ell = int(n/m)                  # m must be a divisor of n
    for k in range(1,m+1):
        r_signal = np.roll(y, k*ell)      
        rot_arr.append(r_signal)
    # the "j" part
    sy = y[:][::-1]                 # Reflecting the signal
    # rotating
    for k in range(1,m+1):
        rjs_signal = np.roll(sy, j+k*ell)  
        rot_arr.append(rjs_signal)
    signal_sym = (1/(2*m))*sum(rot_arr)
    return signal_sym

def is_vector_real(v):
    cond = True
    thd = 1e-13;
    for x in v:
        if np.abs(np.imag(x))>thd:
            cond = False
    return cond

#===========================================================================

# Parameters
n = 15 # make it only odd, so the DC will be in the middle
m = 5 # a divisor of n
j = 2 # the shifted reflection index: can be 0,...,n/m-1

# for generating data set
noise_snr = 5               # desired SNR in dB
num_references = 5          # number of signals
rho = np.random.rand(n);    # the distribution
rho = rho/(sum(rho))

# getting the symmetric signal
signal = make_sym_signal(n, m, j)
plot_discrete(signal,'Signal x')
#signal = make_sym_signal(n, m, j+1)
#plot_discrete(signal,'(real) Signal x2')

# -------- sainity FFT check -----------------
F = np.fft.fft(np.identity(n))
# e = np.linalg.norm(F.dot(signal) - np.fft.fft(signal))
# Fsignal  =  np.fft.fft(signal)
        
#plt.plot(np.abs(Fsignal))
#plt.show()
#plot_discrete(np.abs(Fsignal),'Abs value Fx')
# -------------------------------------------

centered_fft = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(signal)))
#result = (np.fft.fft(np.fft.ifftshift(signal)))
# for x in range(np.size(result)):
#     plt.polar([0,np.angle(result[x])],[0,np.abs(result[x])],marker='o', alpha=0.9)
#     plt.show()
# e = np.linalg.norm(result.imag)

plot_discrete(np.abs(centered_fft),'Abs value of Fx')
# plot_discrete(np.angle(centered_fft),'Angle values of Fx')
print(is_vector_real(centered_fft))
print(np.angle(centered_fft[2]),np.angle(centered_fft[12]))
#is_real = np.isreal(result)
#all_real = np.all(is_real)
#print("All elements are real:", all_real)


# ------------ moments calculations
c_s = circulant(signal)
# print(circ_matrix)

M1 = c_s.dot(rho)
M2 = c_s.dot(np.diag(rho)).dot(np.transpose(c_s))
S = F.dot(M2).dot(np.matrix(F).H) # np.abs(np.fft.fft(M2))
#plt.plot(np.abs(np.fft.fft(np.ones(10))))
#plt.show()

data = generate_data(num_references, signal, noise_snr, rho)
aligned_data = align_data(data)

# Plot original and aligned data
# ipywidgets.interact(plot_data,snr = (0,5,.05))

# plot_data(data, 'Original Data')
# plot_data(aligned_data, 'Aligned Data')

