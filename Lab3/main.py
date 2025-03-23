import numpy as np

def seidel_method(A, b, eps=1e-4, max_iter=100):
    """
    Розв'язок системи A x = b методом Зейделя.
    A - np.array(n,n)
    b - np.array(n)
    eps - точність для умови зупинки
    max_iter - максимальна кількість ітерацій
    Повертає (x, num_iter, history_residuals)
    """
    n = len(A)
    x = np.zeros(n)  # початкове наближення x^(0)
    history_res = []

    for k in range(max_iter):
        x_old = x.copy()
        # послідовне оновлення x_i
        for i in range(n):
            s1 = 0.0
            s2 = 0.0
            # сума по j < i  -> викор. вже "оновлені" x_j
            for j in range(i):
                s1 += A[i,j] * x[j]
            # сума по j > i  -> викор. "старі" x_j
            for j in range(i+1, n):
                s2 += A[i,j] * x_old[j]
            x[i] = (b[i] - s1 - s2)/A[i,i]
        
        # обчислимо нев’язку r^{(k+1)}
        r = b - A.dot(x)
        history_res.append(r)

        # перевірка умови зупинки
        if np.linalg.norm(x - x_old, ord=np.inf) < eps:
            return x, (k+1), history_res

    return x, max_iter, history_res  # якщо за max_iter не досягли eps

# ----------------- Демонстрація ---------------------
if __name__ == "__main__":
    A = np.array([
        [2.12, 0.42, 1.34, 0.88],
        [0.42, 3.95, 1.87, 0.43],
        [1.34, 1.87, 2.98, 0.46],
        [0.88, 0.43, 0.46, 4.44]
    ], dtype=float)
    b = np.array([11.172, 0.115, 0.009, 9.349], dtype=float)

    x_sol, iters, res_hist = seidel_method(A, b, eps=1e-4, max_iter=50)

    print(f"Знайдено x = {x_sol} за {iters} ітерацій.")
    print("Остання нев'язка:", res_hist[-1])
    print("Норма останньої нев'язки:", np.linalg.norm(res_hist[-1], ord=np.inf))
