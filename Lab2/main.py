def gauss_solve(A, b):
    """
    Розв'язує систему A x = b методом Гаусса (без вибору головного елемента).
    A - список списків (матриця n×n),
    b - список (вектор правої частини довжиною n).
    Повертає розв'язок x як список довжиною n.
    """
    import copy
    mat = copy.deepcopy(A)
    vec = copy.deepcopy(b)
    n = len(mat)

    # Прямий хід
    for k in range(n - 1):
        # Перевірка, щоб діагональний елемент не був занадто малим
        if abs(mat[k][k]) < 1.0e-14:
            raise ValueError("Дуже малий (або нульовий) дільник на діагоналі!")
        for i in range(k + 1, n):
            c = mat[i][k] / mat[k][k]
            for j in range(k, n):
                mat[i][j] -= c * mat[k][j]
            vec[i] -= c * vec[k]

    # Зворотний хід
    x = [0.0]*n
    for i in reversed(range(n)):
        s = vec[i]
        for j in range(i+1, n):
            s -= mat[i][j] * x[j]
        if abs(mat[i][i]) < 1.0e-14:
            raise ValueError("Нульовий дільник на діагоналі під час зворотного ходу!")
        x[i] = s / mat[i][i]
    return x

# Приклад використання:
if __name__ == "__main__":
    A = [
        [2.12, 0.42, 1.34, 0.88],
        [0.42, 3.95, 1.87, 0.43],
        [1.34, 1.87, 2.98, 0.46],
        [0.88, 0.43, 0.46, 4.44]
    ]
    b = [11.172, 0.115, 0.009, 9.349]

    sol = gauss_solve(A, b)
    print("Розв'язок (x1, x2, x3, x4):", sol)

    # Перевіримо нев'язку
    import numpy as np
    A_np = np.array(A, dtype=float)
    b_np = np.array(b, dtype=float)
    x_np = np.array(sol, dtype=float)
    resid = b_np - A_np.dot(x_np)
    print("Нев'язка:", resid)
