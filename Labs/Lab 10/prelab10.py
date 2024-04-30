import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre

def eval_legendre(n, x):
    p = [1]  
    if n > 0:
        p.append(x)  
    
    for i in range(2, n + 1):
        p.append(((2 * i - 1) * x * p[i - 1] - (i - 1) * p[i - 2]) / i)
    
    return p

n = 5
x_values = np.linspace(-1, 1, 100)
legendre_values = [eval_legendre(n, x) for x in x_values]

plt.figure(figsize=(10, 6))
for i in range(n + 1):
    scipy_legendre_values = legendre(i)(x_values)
    plt.plot(x_values, scipy_legendre_values, linestyle='--', label=f'SciPy: n={i}')

plt.title('Legendre Polynomials Comparison')
plt.xlabel('x')
plt.ylabel('Legendre Polynomials')
plt.legend()
plt.grid(True)
plt.show()
