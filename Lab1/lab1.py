import math

def f(x):
    return -x**4 + 3*x**3 - 2*x + 4

def f_prime(x):
    """Похідна f(x)."""
    return -4*x**3 + 9*x**2 - 2

# 1) Відокремлення коренів (сканування на проміжку від -5 до 5 кроком 1.0)
def find_intervals(a, b, step=1.0):
    intervals = []
    x_left = a
    f_left = f(x_left)
    while x_left < b:
        x_right = x_left + step
        f_right = f(x_right)
        if f_left * f_right < 0:
            intervals.append((x_left, x_right))
        x_left = x_right
        f_left = f_right
    return intervals

# 2) Метод бісекції
def bisection(a, b, eps=1e-5):
    fa = f(a)
    fb = f(b)
    if fa*fb > 0:
        raise ValueError("На відрізку немає зміни знака.")
    iter_num = 0
    while (b - a)/2 > eps:
        iter_num += 1
        m = (a + b)/2.0
        fm = f(m)
        print(f"Бісекція, ітерація {iter_num}: a = {a:.6f}, b = {b:.6f}, m = {m:.6f}, f(m) = {fm:.6e}, (b-a)/2 = {(b-a)/2:.6e}")
        if abs(fm) < eps:
            return m
        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
    return (a + b)/2.0

# 3) Метод хорд
def secant(a, b, eps=1e-5):
    fa = f(a)
    fb = f(b)
    if fa*fb > 0:
        raise ValueError("На відрізку немає зміни знака.")
    x0, x1 = a, b
    f0, f1 = fa, fb
    iter_num = 0
    while abs(x1 - x0) > eps:
        iter_num += 1
        x2 = x1 - f1*(x1 - x0)/(f1 - f0)
        f2 = f(x2)
        print(f"Хорди, ітерація {iter_num}: x0 = {x0:.6f}, x1 = {x1:.6f}, x2 = {x2:.6f}, f(x2) = {f2:.6e}, |x2-x1| = {abs(x2-x1):.6e}")
        if abs(x2 - x1) < eps:
            return x2
        x0, x1 = x1, x2
        f0, f1 = f1, f2
    return x1

# 4) Метод Ньютона
def newton(x0, eps=1e-5, max_iter=100):
    iter_num = 0
    while iter_num < max_iter:
        iter_num += 1
        fx = f(x0)
        fpx = f_prime(x0)
        if abs(fpx) < 1e-14:
            raise ZeroDivisionError("Похідна близька до нуля.")
        x1 = x0 - fx/fpx
        print(f"Ньютон, ітерація {iter_num}: x = {x0:.6f}, f(x) = {fx:.6e}, f'(x) = {fpx:.6e}, x_new = {x1:.6f}, |x_new-x| = {abs(x1-x0):.6e}")
        if abs(x1 - x0) < eps:
            return x1
        x0 = x1
    return x0


intervals = find_intervals(-5, 5, step=1.0)
print("Знайдені проміжки зі зміною знака:")
for (a_int, b_int) in intervals:
    print(f"[{a_int}, {b_int}]  f(a)={f(a_int):.3f}, f(b)={f(b_int):.3f}")

eps = 1e-5
for (a_int, b_int) in intervals:
    print(f"\n=== Уточнення на відрізку [{a_int}, {b_int}] ===")
    try:
        root_bis = bisection(a_int, b_int, eps)
        root_sec = secant(a_int, b_int, eps)
        root_newton = newton((a_int + b_int)/2, eps)
        
        print(f"Метод бісекції: x = {root_bis:.6f}, f(x) = {f(root_bis):.6e}")
        print(f"Метод хорд:    x = {root_sec:.6f}, f(x) = {f(root_sec):.6e}")
        print(f"Метод Ньютона: x = {root_newton:.6f}, f(x) = {f(root_newton):.6e}")
    except Exception as e:
        print(f"Помилка: {e}")
