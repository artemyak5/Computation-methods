import numpy as np
import matplotlib.pyplot as plt

# Параметри
a, b = 0, 4
n = 10
x_nodes = np.linspace(a, b, n+1)
y_nodes = x_nodes * np.sqrt(x_nodes)

# Функція f(x)
def f(x):
    return x * np.sqrt(x)

# 1. Розділені різниці (коефіцієнти полінома Ньютона)
def divided_differences(x, y):
    n = len(x)
    dd = np.zeros((n, n))
    dd[:, 0] = y
    for j in range(1, n):
        for i in range(n - j):
            dd[i, j] = (dd[i+1, j-1] - dd[i, j-1]) / (x[i+j] - x[i])
    return dd[0]

coeffs = divided_differences(x_nodes, y_nodes)

# 2. Оцінка полінома Ньютона
def newton_poly(x, x_nodes, coeffs):
    n = len(coeffs)
    result = coeffs[-1]
    for k in range(n-2, -1, -1):
        result = result * (x - x_nodes[k]) + coeffs[k]
    return result

# 3. Кубічний сплайн (натуральний)
h = np.diff(x_nodes)
A = np.zeros((n+1, n+1))
b_vec = np.zeros(n+1)
A[0,0] = 1
A[-1,-1] = 1
for i in range(1, n):
    A[i, i-1] = h[i-1]
    A[i, i]   = 2*(h[i-1] + h[i])
    A[i, i+1] = h[i]
    b_vec[i]  = 6*((y_nodes[i+1] - y_nodes[i]) / h[i] - (y_nodes[i] - y_nodes[i-1]) / h[i-1])

m = np.linalg.solve(A, b_vec)

def spline(x):
    i = np.searchsorted(x_nodes, x) - 1
    i = np.clip(i, 0, n-1)
    xi, xi1 = x_nodes[i], x_nodes[i+1]
    hi = xi1 - xi
    Mi, Mi1 = m[i], m[i+1]
    yi, yi1 = y_nodes[i], y_nodes[i+1]
    return (Mi*(xi1-x)**3 + Mi1*(x-xi)**3)/(6*hi) + (yi - Mi*hi**2/6)*(xi1-x)/hi + (yi1 - Mi1*hi**2/6)*(x-xi)/hi

# 4. Густий сітковий розрахунок
x_dense = np.linspace(a, b, 500)
y_true = f(x_dense)
y_poly = newton_poly(x_dense, x_nodes, coeffs)
y_spline = np.array([spline(x) for x in x_dense])
error_poly = np.abs(y_poly - y_true)

plt.figure()
plt.plot(x_dense, y_true, label='f(x)')
plt.plot(x_dense, y_poly, label='Newton polynomial')
plt.scatter(x_nodes, y_nodes, marker='o')
plt.title('Polynomial interpolation (Newton)')
plt.legend()

plt.figure()
plt.plot(x_dense, y_true, label='f(x)')
plt.plot(x_dense, y_spline, label='Cubic spline')
plt.scatter(x_nodes, y_nodes, marker='o')
plt.title('Cubic spline interpolation')
plt.legend()

plt.figure()
plt.plot(x_dense, error_poly)
plt.title('Error |P_n(x)-f(x)|')
plt.tight_layout()

plt.show()
