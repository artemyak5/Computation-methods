import math

def cholesky_decomposition(A):
    n = len(A)
    # Ініціалізація T як нульової нижньої трикутної матриці
    T = [[0.0]*n for _ in range(n)]
    for i in range(n):
        # 1) t_{ii}
        s = A[i][i]
        for k in range(i):
            s -= T[i][k]**2
        if s <= 0:
            raise ValueError("Матриця не є додатньо визначеною (від’ємний діагональний елемент).")
        T[i][i] = math.sqrt(s)
        # 2) t_{ji} (j>i)
        for j in range(i+1, n):
            s = A[j][i]
            for k in range(i):
                s -= T[j][k]*T[i][k]
            T[j][i] = s / T[i][i]
    return T

def forward_substitution(T, b):
    """Розв'язує T*y = b, де T - нижня трикутна матриця."""
    n = len(T)
    y = [0.0]*n
    for i in range(n):
        s = b[i]
        for k in range(i):
            s -= T[i][k]*y[k]
        y[i] = s / T[i][i]
    return y

def backward_substitution(Tt, y):
    """Розв'язує T^T*x = y (Tt = T^T - верхня трикутна)."""
    n = len(Tt)
    x = [0.0]*n
    for i in reversed(range(n)):
        s = y[i]
        for k in range(i+1, n):
            s -= Tt[i][k]*x[k]
        x[i] = s / Tt[i][i]
    return x

# Вхідні дані (матриця A та вектор b):
A = [
    [2.12, 0.42, 1.34, 0.88],
    [0.42, 3.95, 1.87, 0.43],
    [1.34, 1.87, 2.98, 0.46],
    [0.88, 0.43, 0.46, 4.44]
]
b = [11.172, 0.115, 0.009, 9.349]

# 1) Отримуємо матрицю T (нижню трикутну)
try:
    T = cholesky_decomposition(A)
    print("Матриця T (нижня трикутна):")
    for row in T:
        print(row)
    
    # 2) Розв'язуємо T*y = b
    y = forward_substitution(T, b)
    print("\nВектор y =")
    print(y)
    
    # 3) Розв'язуємо T^T*x = y
    # Побудуємо T^T (верхня трикутна)
    Tt = [[T[j][i] for j in range(len(T))] for i in range(len(T))]
    x = backward_substitution(Tt, y)
    
    print("\nРозв'язок системи (x):")
    for i, xi in enumerate(x):
        print(f"x[{i+1}] = {xi:.6f}")
    
    # Перевірка нев’язки r = b - A*x
    r = []
    for i in range(len(A)):
        ax_i = sum(A[i][j]*x[j] for j in range(len(x)))
        r.append(b[i] - ax_i)
    print("\nВектор нев'язки r = b - A*x:")
    for i, ri in enumerate(r):
        print(f"r[{i+1}] = {ri:.6e}")

except ValueError as e:
    print("Помилка:", e)