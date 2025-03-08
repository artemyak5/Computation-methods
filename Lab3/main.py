import numpy as np

def simple_iteration(A, b, tol=1e-6, max_iter=200):
    n = len(b)
    x = np.zeros(n)  # початкове наближення (0,0,0,0)
    iteration_history = []  # для виведення кожної ітерації

    for k in range(max_iter):
        x_new = np.zeros(n)
        # Обчислення нового вектора за формулами:
        x_new[0] = (b[0] - (A[0,1]*x[1] + A[0,2]*x[2] + A[0,3]*x[3])) / A[0,0]
        x_new[1] = (b[1] - (A[1,0]*x[0] + A[1,2]*x[2] + A[1,3]*x[3])) / A[1,1]
        x_new[2] = (b[2] - (A[2,0]*x[0] + A[2,1]*x[1] + A[2,3]*x[3])) / A[2,2]
        x_new[3] = (b[3] - (A[3,0]*x[0] + A[3,1]*x[1] + A[3,2]*x[2])) / A[3,3]

        iteration_history.append((k+1, x_new.copy(), np.linalg.norm(x_new - x, ord=np.inf)))
        if np.linalg.norm(x_new - x, ord=np.inf) < tol:
            return x_new, iteration_history
        x = x_new.copy()
    return x, iteration_history

# Задання системи (варіант 4)
A = np.array([
    [2.12,  0.42,  1.34,  0.88],
    [0.42,  3.95,  1.87,  0.43],
    [1.34,  1.87,  2.98,  0.46],
    [0.88,  0.43,  0.46,  4.44]
])
b = np.array([11.172, 0.115, 0.009, 9.349])

solution, history = simple_iteration(A, b, tol=1e-6, max_iter=200)

# Вивід результатів і кожної ітерації
print("Отриманий розв'язок:")
print(solution)
print("\nІтерації (номер, x^(k+1), ||x^(k+1)-x^(k)||∞):")
for it_num, x_val, diff in history:
    print(f"Ітерація {it_num}: x = {x_val}, критерій = {diff}")

# Обчислення нев'язки r = b - A*x
r = b - A.dot(solution)
print("\nВектор нев'язки r =")
print(r)
print("||r||∞ =", np.linalg.norm(r, ord=np.inf))
