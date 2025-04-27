def seidel_method_no_lib(A, b, eps=1e-4, max_iter=100):
    n = len(A)
    # Початкове наближення x^(0): усе нулі
    x = [0.0]*n
    
    # Список для збереження векторів нев'язки
    history_res = []
    
    # Функція для множення матриці A на вектор x: r_i = Σ_j A[i][j]*x[j]
    def mat_vec_mult(mat, vec):
        m = len(mat)
        res = [0.0]*m
        for i in range(m):
            s = 0.0
            for j in range(len(vec)):
                s += mat[i][j]*vec[j]
            res[i] = s
        return res
    
    # Функція для обчислення вектора нев'язки r = b - A*x
    def residual(A, x, b):
        Ax = mat_vec_mult(A, x)
        r = [b[i] - Ax[i] for i in range(len(b))]
        return r
    
    # Функція для обчислення "норми" різниці (тут використовуємо максимум за модулем)
    def inf_norm_diff(vec1, vec2):
        mx = 0.0
        for i in range(len(vec1)):
            diff = abs(vec1[i] - vec2[i])
            if diff > mx:
                mx = diff
        return mx
    
    # Визначимо r^(0)
    r0 = residual(A, x, b)
    history_res.append(r0)
    
    # Ітераційний процес
    for k in range(max_iter):
        x_old = x[:]
        
        # Обчислюємо x_i^{(k+1)} за формулами методу Зейделя
        for i in range(n):
            # sum1: сума A[i][j]*x[j] для j < i (вже оновлені x_j)
            sum1 = 0.0
            for j in range(i):
                sum1 += A[i][j]*x[j]
            # sum2: сума A[i][j]*x_old[j] для j > i (старі x_j)
            sum2 = 0.0
            for j in range(i+1, n):
                sum2 += A[i][j]*x_old[j]
            
            x[i] = (b[i] - sum1 - sum2)/A[i][i]
        
        # Обчислимо вектор нев'язки після (k+1)-ї ітерації
        r_kplus1 = residual(A, x, b)
        history_res.append(r_kplus1)
        
        # Перевірка зупинки за зміною x
        if inf_norm_diff(x, x_old) < eps:
            return (x, k+1, history_res)
    
    # Якщо не вклалися у max_iter
    return (x, max_iter, history_res)


if __name__ == "__main__":
    # Припустимо, ми маємо 4×4-систему (варіант 4)
    A = [
        [2.12, 0.42, 1.34, 0.88],
        [0.42, 3.95, 1.87, 0.43],
        [1.34, 1.87, 2.98, 0.46],
        [0.88, 0.43, 0.46, 4.44]
    ]
    b = [11.172, 0.115, 0.009, 9.349]
    
    # Виклик методу Зейделя
    x_sol, iters, res_history = seidel_method_no_lib(A, b, eps=1e-4, max_iter=50)
    
    print(f"Отриманий розв'язок x після {iters} ітерацій:")
    for i, val in enumerate(x_sol, start=1):
        print(f"x_{i} = {val:.5f}")
    
    print("\nНев'язка на останній ітерації:")
    r_last = res_history[-1]
    for i, val in enumerate(r_last, start=1):
        print(f"r_{i} = {val:.5e}")
    
    # for k, rv in enumerate(res_history):
    #     print(f"Ітерація {k}: r = {rv}")
