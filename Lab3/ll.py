def simple_iteration_method(A, b, eps=1e-4, max_iter=1000):
    n = len(A)
    x_old = [0.0] * n
    for k in range(max_iter):
        x_new = [0.0] * n
        for i in range(n):
            sigma = 0.0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x_old[j]
            x_new[i] = (b[i] - sigma) / A[i][i]
        diff = max(abs(x_new[i] - x_old[i]) for i in range(n))
        r = [b[i] - sum(A[i][j]*x_new[j] for j in range(n)) for i in range(n)]
        print(f"Iter {k+1}: x = {[round(val,5) for val in x_new]}, ||dx|| = {diff:.5e}, residual = {[round(abs(ri),5) for ri in r]}")
        if diff < eps:
            break
        x_old = x_new
if __name__ == "__main__":
    A = [
        [2.12, 0.42, 1.34, 0.88],
        [0.42, 3.95, 1.87, 0.43],
        [1.34, 1.87, 2.98, 0.46],
        [0.88, 0.43, 0.46, 4.44]
    ]
    b = [11.172, 0.115, 0.009, 9.349]
    simple_iteration_method(A, b, eps=1e-4)
