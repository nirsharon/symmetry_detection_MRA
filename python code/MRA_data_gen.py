# generating plain MRA data (cyclic group)

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import circulant
from EM import *
#import ipywidgets

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
    n_order =  np.size(signal)
    # x = np.linspace(0, 2 * np.pi, n_order) #sample uniform [0,2pi]
    y = signal #np.sin(x)
    values = range(0,n_order)
    group_sampls = random.choices(values, rho, k=num_references);
    
    for el in range(num_references):
        shift_amount = group_sampls[el] # np.random.randint(0,n_order-1)
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

 # Plot data varies snr <======== ignore for now
 # def plot_data(snr)
 #    num_references1 = 5
 #    n_order1 = 100
 #    title = 'varied SNR'
 #    data = generate_data(num_references1, n_order1, snr)
 #    for y in data:
 #        plt.plot(y, alpha=0.5)
 #    plt.title(title)
 #    plt.show()

def plot_discrete(vec, nameit):
    
    n = np.size(vec)
    plt.stem(range(n), vec, use_line_collection=True)
    #plt.xlabel('n')
    #plt.xticks(np.arange(LL, UL, 1))
    #plt.yticks([0, 1])
    plt.ylabel('value')
    plt.title(nameit)
    plt.show()


# Parameters
num_references = 5
n_order = 10
m_sym = 5 # a divisor of n_order
noise_var = 0.1
noise_snr = 5  # Desired SNR in dB

rho = np.random.rand(n_order);
rho = rho/(sum(rho))

# Generate and align data
signal = make_signal(n_order, m_sym);

plt.plot(signal, alpha=0.9)
plt.show()

plot_discrete(signal,'Signal x')


c_s = circulant(signal)
# print(circ_matrix)

M1 = c_s.dot(rho)
M2 = c_s.dot(np.diag(rho)).dot(np.transpose(c_s))


F = np.fft.fft(np.identity(n_order))
e = np.linalg.norm(F.dot(signal) - np.fft.fft(signal))
Fsignal   =    np.fft.fft(signal);
M1_hat = np.fft.fft(M1);
        
plt.plot(np.abs(Fsignal))
plt.show()
plot_discrete(np.abs(Fsignal),'Abs value Fx')
plot_discrete(np.abs(M1_hat),'First moment')


S = F.dot(M2).dot(np.matrix(F).H) # np.abs(np.fft.fft(M2))
plt.imshow(np.abs(S))
power_spectrum = np.diagonal(S)
plot_discrete(np.abs(Fsignal),'Abs value Fx')
#plt.plot(np.abs(np.fft.fft(np.ones(10))))
#plt.show()

data = generate_data(num_references, signal, 50, rho)
aligned_data = align_data(data)

# Plot original and aligned data
# ipywidgets.interact(plot_data,snr = (0,5,.05))

# plot_data(data, 'Original Data')
# plot_data(aligned_data, 'Aligned Data')

MRA_data = generate_data(num_references, signal, noise_snr, rho)
em_algorithm(MRA_data)