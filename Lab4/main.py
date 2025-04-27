import math

# Матриця варіанта 4
A = [
    [6.30, 1.07, 0.99, 1.20],
    [1.07, 4.12, 1.30, 0.16],
    [0.99, 1.30, 5.48, 2.10],
    [1.20, 0.16, 2.10, 6.06]
]

y0 = [1.0, 0.0, 0.0, 0.0]

def mat_vec_mul(M, v):
    n = len(v)
    return [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]

y1 = mat_vec_mul(A, y0)
y2 = mat_vec_mul(A, y1)
y3 = mat_vec_mul(A, y2)
y4 = mat_vec_mul(A, y3)

C = [[y3[i], y2[i], y1[i], y0[i]] for i in range(4)]
b = [-y4[i] for i in range(4)]

def gaussian_elimination(aug):
    n = len(aug)
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(aug[r][i]))
        aug[i], aug[max_row] = aug[max_row], aug[i]
        pivot = aug[i][i]
        for j in range(i+1, n):
            factor = aug[j][i] / pivot
            for k in range(i, n+1):
                aug[j][k] -= factor * aug[i][k]
    x = [0] * n
    for i in range(n-1, -1, -1):
        s = aug[i][n] - sum(aug[i][j] * x[j] for j in range(i+1, n))
        x[i] = s / aug[i][i]
    return x

aug_matrix = [C[i] + [b[i]] for i in range(4)]
p = gaussian_elimination(aug_matrix)

print("Коефіцієнти p_i (метод Крилова):")
for idx, pi in enumerate(p, 1):
    print(f"p{idx} = {pi:.6f}")

def poly(x):
    return x**4 + p[0]*x**3 + p[1]*x**2 + p[2]*x + p[3]

def poly_der(x):
    return 4*x**3 + 3*p[0]*x**2 + 2*p[1]*x + p[2]

def newton(x0, iterations=20):
    x = x0
    for _ in range(iterations):
        x = x - poly(x)/poly_der(x)
    return x

initial_guesses = [2.5, 4.0, 6.0, 9.0]
roots = [newton(g) for g in initial_guesses]

print("\nНаближені власні числа (метод Ньютона):")
for r in roots:
    print(f"{r:.6f}")
