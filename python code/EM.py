import numpy as np
from scipy.ndimage import shift
import matplotlib.pyplot as plt

def initialize_parameters(signals, length):
    n, d = signals.shape
    mean_signal = np.mean(signals, axis=0)
    variance = np.var(signals)
    return mean_signal, variance

def e_step(signals, mean_signal, variance):
    n, d = signals.shape
    responsibilities = np.zeros((n, d))
    for i in range(n):
        for shift_amount in range(d):
            shifted_signal = np.roll(signals[i], shift_amount)
            responsibilities[i, shift_amount] = np.exp(-np.sum((shifted_signal - mean_signal)**2) / (2 * variance))
    responsibilities /= responsibilities.sum(axis=1, keepdims=True)
    return responsibilities

def m_step(signals, responsibilities):
    n, d = signals.shape
    mean_signal = np.zeros(d)
    variance = 0
    for i in range(n):
        for shift_amount in range(d):
            shifted_signal = np.roll(signals[i], shift_amount)
            mean_signal += responsibilities[i, shift_amount] * shifted_signal
    mean_signal /= n
    for i in range(n):
        for shift_amount in range(d):
            shifted_signal = np.roll(signals[i], shift_amount)
            variance += responsibilities[i, shift_amount] * np.sum((shifted_signal - mean_signal)**2)
    variance /= (n * d)
    return mean_signal, variance

def em_algorithm(signals, max_iter=100, tol=1e-4):
    mean_signal, variance = initialize_parameters(signals, signals.shape[1])
    log_likelihoods = []
    for i in range(max_iter):
        responsibilities = e_step(signals, mean_signal, variance)
        mean_signal, variance = m_step(signals, responsibilities)
        log_likelihood = np.sum(np.log(np.sum([np.exp(-np.sum((np.roll(signals[i], shift_amount) - mean_signal)**2) / (2 * variance)) for shift_amount in range(signals.shape[1])], axis=0)))
        log_likelihoods.append(log_likelihood)
        if len(log_likelihoods) > 1 and np.abs(log_likelihoods[-1] - log_likelihoods[-2]) < tol:
            break
    return mean_signal, variance, log_likelihoods

# Example usage
signals = np.random.randn(100, 10)  # Replace with your signal data
mean_signal, variance, log_likelihoods = em_algorithm(signals)
print("Recovered Signal:", mean_signal)
print("Variance:", variance)
plt.plot(mean_signal)
plt.show()