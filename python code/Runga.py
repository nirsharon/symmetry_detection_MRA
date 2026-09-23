import numpy as np
import matplotlib.pyplot as plt

# Function to compute barycentric interpolation
def barycentric_interpolation(nodes, values, x):
    n = len(nodes)
    weights = np.ones(n)
    for j in range(n):
        for k in range(n):
            if j != k:
                weights[j] /= (nodes[j] - nodes[k])
    
    numerator = np.zeros_like(x)
    denominator = np.zeros_like(x)
    exact = np.zeros_like(x, dtype=bool)
    
    for j in range(n):
        diff = x - nodes[j]
        mask = diff == 0
        exact[mask] = True
        numerator += weights[j] * values[j] / diff
        denominator += weights[j] / diff
    
    result = numerator / denominator
    result[exact] = values[np.argmax(exact)]
    return result

# Function
f = lambda x: 1 / (x**2 + 1)
# f = lambda x: np.cos(np.pi / 5 * x)

MAXDEG = 25  # Maximal degree of the interpolating polynomial
a, b = -5, 5  # Interpolation interval

xx = np.linspace(a, b, 1000)  # Plotting grid

for n in range(1, MAXDEG + 1):
    # Plot the function
    plt.clf()
    plt.plot(xx, f(xx), 'b', label="Original Function")
    
    # Interpolate in equidistant nodes
    nodes = np.linspace(a, b, n + 1)
    # Chebyshev nodes (optional)
    #nodes = (a + b) / 2 + ((b - a) / 2) * np.cos((2 * np.arange(1, n + 1) - 1) / (2 * n) * np.pi)
    intf = barycentric_interpolation(nodes, f(nodes), xx)
    
    # Plot the interpolation
    plt.plot(nodes, f(nodes), 'r.', markersize=15, label="Interpolation Nodes")
    plt.plot(xx, intf, 'r', label="Interpolated Polynomial")
    plt.ylim([-1, 2])
    
    # Add error as title
    err = np.max(np.abs(intf - f(xx)))
    plt.title(f"Polynomial interpolation degree = {n}, err = {err:.3e}")
    plt.legend()
    plt.pause(0.5)