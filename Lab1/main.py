import math

def f(x):
    return -x**4 + 3*x**3 - 2*x + 4

def df(x):
    # Похідна від f(x)
    return -4*x**3 + 9*x**2 - 2

# 1. Метод бісекції
def bisection(a, b, eps=1e-6, max_iter=100):
    print("\n=== Метод бісекції ===")
    for i in range(max_iter):
        c = (a + b)/2
        fc = f(c)
        print(f"Ітерація {i}: c = {c:.6f}, f(c) = {fc:.6e}, інтервал = [{a:.6f}, {b:.6f}]")
        if abs(fc) < eps:
            return c, i+1
        if f(a)*fc < 0:
            b = c
        else:
            a = c
        if (b - a) < eps:
            return (a + b)/2, i+1
    return (a + b)/2, max_iter

# 2. Метод хорд (січних)
def chord(a, b, eps=1e-6, max_iter=100):
    print("\n=== Метод хорд ===")
    x_old = a
    for i in range(max_iter):
        fa = f(a)
        fb = f(b)
        x_new = (a*fb - b*fa)/(fb - fa)  # Формула хорди
        fx_new = f(x_new)
        print(f"Ітерація {i}: x_new = {x_new:.6f}, f(x_new) = {fx_new:.6e}")
        if abs(fx_new) < eps:
            return x_new, i+1
        # Визначаємо, куди «зсунути» інтервал
        if fa * fx_new < 0:
            b = x_new
        else:
            a = x_new
        if abs(x_new - x_old) < eps:
            return x_new, i+1
        x_old = x_new
    return x_new, max_iter

# 3. Метод Ньютона (дотичних)
def newton(x0, eps=1e-6, max_iter=100):
    print("\n=== Метод Ньютона ===")
    for i in range(max_iter):
        fx0 = f(x0)
        dfx0 = df(x0)
        if abs(dfx0) < 1e-12:
            print("Похідна надто мала, можливе зупинення.")
            return x0, i
        x1 = x0 - fx0/dfx0
        print(f"Ітерація {i}: x = {x1:.6f}, f(x) = {f(x1):.6e}")
        if abs(x1 - x0) < eps:
            return x1, i+1
        x0 = x1
    return x0, max_iter

# Основний виклик для відрізків:
intervals = [(-2, -1), (2, 3)]
for (left, right) in intervals:
    print(f"\n=== УТОЧНЕННЯ КОРЕНЯ на [{left}, {right}] ===")
    # Бісекція
    root_bis, it_bis = bisection(left, right)
    print(f"Корінь (бісекція) ≈ {root_bis:.6f}, ітерацій: {it_bis}")
    # Хорди
    root_ch, it_ch = chord(left, right)
    print(f"Корінь (хорди)    ≈ {root_ch:.6f}, ітерацій: {it_ch}")
    # Ньютона: початкове наближення візьмемо середину
    x0 = (left + right)/2
    root_newt, it_newt = newton(x0)
    print(f"Корінь (Ньютона)  ≈ {root_newt:.6f}, ітерацій: {it_newt}")
    print("-"*50)
