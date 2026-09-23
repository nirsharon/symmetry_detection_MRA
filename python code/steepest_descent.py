import numpy as np
import matplotlib.pyplot as plt

# Matrix 1
A = np.array([[1, 0],
              [0, 0.1]])
# Matrix 2
A = np.array([[1, 0],
              [0, 0.8]])

# Matrix 3
#A = np.random.randn(2, 2)
#A = A @ A.T
eigenvalues = np.linalg.eigvals(A)
print(f'The eigenvalues are {eigenvalues[0]:.4f} and {eigenvalues[1]:.4f}')
b = np.array([0.1, 0.1])

def f(x):
    return (x.T @ A @ x) / 2 - b.T @ x

T = 100
N = 100
t = np.linspace(-T, T, N)

X, Y = np.meshgrid(t, t)
Z = np.zeros_like(X)

for j in range(Z.size):
    Z.flat[j] = f(np.array([X.flat[j], Y.flat[j]]))

fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, facecolors=plt.cm.jet(Z / Z.max()), edgecolor='none')
plt.show()

fig = plt.figure(2)
plt.clf()
plt.contour(X, Y, Z)
plt.axis('equal')
plt.axis([-T, T, -T, T])

# fig = plt.figure(2)
# plt.hold(True)
# x0 = np.array([T/10, T])  # Slow
# x0 = np.array([0, T])  # Super fast
x0 = np.array([T, T])  # Fast
xi = x0
ri = 100
iter = 0
plt.scatter(xi[0], xi[1], 5, 'black')
plt.text(xi[0], xi[1], f'{iter}')
print('Press any key to start')
input()
while np.abs(ri).max() > 1.0e-3:
    ri = b - A @ xi
    alphai = (ri.T @ ri) / (ri.T @ A @ ri)
    np_ = xi + alphai * ri
    print(f'iter={iter}\t xi=[{xi[0]:+.3e},{xi[1]:+.3e}] \t ri=[{ri[0]:+.1e},{ri[1]:+.1e}] \t ||ri||={np.linalg.norm(ri):.1e}')
    iter += 1
    plt.scatter(np_[0], np_[1], 5, 'black')
    plt.plot([xi[0], np_[0]], [xi[1], np_[1]])
    plt.text(np_[0], np_[1], f'{iter}')
    xi = np_
    plt.pause(1)
#plt.hold(False)
print(f'Final xi=[{xi[0]:.4f}, {xi[1]:.4f}]')
plt.show()