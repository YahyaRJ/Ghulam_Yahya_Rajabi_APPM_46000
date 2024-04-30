import numpy as np
import scipy.linalg as scila
import time
import matplotlib.pyplot as plt

def solve_with_normal_eq(A, b):
    x = np.linalg.inv(A.T @ A) @ A.T @ b
    return x

def solve_with_qr(A, b):
    Q, R = np.linalg.qr(A)
    x = np.linalg.solve(R, Q.T @ b)
    return x

def create_rect(N, M):
    ''' this subroutine creates an ill-conditioned rectangular matrix'''
    a = np.linspace(1, 10, M)
    d = 10 ** (-a)

    D2 = np.zeros((N, M))
    for j in range(0, M):
        D2[j, j] = d[j]

    A = np.random.rand(N, N)
    Q1, R = np.linalg.qr(A)
    A = np.random.rand(M, M)
    Q2, R = np.linalg.qr(A)

    B = np.matmul(Q1, D2)
    B = np.matmul(B, Q2)
    return B

def measure_time(N_values, num_rhs, technique='normal_eq'):
    times = []

    for N in N_values:
        A = create_rect(N, num_rhs)  # Use rectangular matrices
        b = np.random.rand(N, num_rhs)

        start = time.time()
        if technique == 'normal_eq':
            x = solve_with_normal_eq(A, b)
        elif technique == 'qr':
            x = solve_with_qr(A, b)
        else:
            raise ValueError("Invalid technique")
        total_time = time.time() - start

        times.append(total_time)

    return times

def driver():
    N_values = [100, 500, 1000, 2000, 4000, 5000]
    num_rhs = 10  # Number of right hand sides

    normal_eq_times = measure_time(N_values, num_rhs, technique='normal_eq')
    qr_times = measure_time(N_values, num_rhs, technique='qr')

    plt.plot(N_values, normal_eq_times, label='Normal Equation')
    plt.plot(N_values, qr_times, label='QR Factorization')
    plt.xlabel('Matrix Size (N)')
    plt.ylabel('Time (seconds)')
    plt.title('Solve Time Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    driver()
