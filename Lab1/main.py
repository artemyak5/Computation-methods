import math

def f(x):
    return -x**4 + 3*x**3 - 2*x + 4

def f_prime(x):
    """Похідна f(x)."""
    return -4*x**3 + 9*x**2 - 2

# 1) Відокремлення коренів (сканування на проміжку, наприклад, від -5 до 5 кроком 0.5)
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
    
    while (b - a)/2 > eps:
        m = (a + b)/2.0
        fm = f(m)
        if abs(fm) < eps:
            return m
        if fa*fm < 0:
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
    while abs(x1 - x0) > eps:
        # Формула хорди
        x2 = x1 - f1*(x1 - x0)/(f1 - f0)
        f2 = f(x2)
        
        # Перевірка збіжності
        if abs(x2 - x1) < eps:
            return x2
        
        # Оновлюємо точки
        x0, x1 = x1, x2
        f0, f1 = f1, f2
    
    return x1

# 4) Метод Ньютона
def newton(x0, eps=1e-5, max_iter=100):
    for i in range(max_iter):
        fx = f(x0)
        fpx = f_prime(x0)
        if abs(fpx) < 1e-14:
            # Уникаємо ділення на нуль
            raise ZeroDivisionError("Похідна близька до нуля.")
        x1 = x0 - fx/fpx
        if abs(x1 - x0) < eps:
            return x1
        x0 = x1
    return x0

# --- ГОЛОВНА ЧАСТИНА ---

# Знайдемо інтервали зі зміною знака
intervals = find_intervals(-5, 5, step=1.0)
print("Знайдені проміжки зі зміною знака:")
for (a_int, b_int) in intervals:
    print(f"[{a_int}, {b_int}]  f(a)={f(a_int):.3f}, f(b)={f(b_int):.3f}")

# Очікуємо, що серед них будуть [-2, -1] і [2, 3]
# Уточнюємо корені на кожному з них
eps = 1e-5
for (a_int, b_int) in intervals:
    print(f"\n=== УТОЧНЕННЯ НА ВІДРІЗКУ [{a_int}, {b_int}] ===")
    try:
        root_bis = bisection(a_int, b_int, eps)
        root_sec = secant(a_int, b_int, eps)
        # Для Ньютона обираємо початок в середині відрізка:
        root_newton = newton((a_int+b_int)/2, eps)
        
        print(f"Метод бісекції: x = {root_bis:.6f}, f(x)={f(root_bis):.6e}")
        print(f"Метод хорд:    x = {root_sec:.6f}, f(x)={f(root_sec):.6e}")
        print(f"Метод Ньютона: x = {root_newton:.6f}, f(x)={f(root_newton):.6e}")
    except Exception as e:
        print(f"Помилка: {e}")
