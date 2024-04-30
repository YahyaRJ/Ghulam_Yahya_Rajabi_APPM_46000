def composite_trapezoidal_rule(a, b, f, N):
    h = (b - a) / N
    result = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        result += f(a + i * h)
    result *= h
    return result

def f(x):
    return x**2

a = 0
b = 1

N = 100

integral = composite_trapezoidal_rule(a, b, f, N)
print(integral)
