import numpy as np

def calculate_derivatives():
    # Define the function
    func = lambda x: np.arcsin(x)

    # Constants
    x_val = np.pi / 4
    actual_derivative = 1 / np.sqrt(1 - x_val**2)

    # Generate h values
    h_values = 0.03 * 2 ** (-1 * np.arange(0, 10, dtype=float))

    # Forward difference method
    print("FORWARD DIFFERENCE METHOD")
    forward_diff_results = calculate_forward_difference(func, h_values, x_val)
    print("Approximations:")
    print(forward_diff_results)
    alpha_val = 1
    print("\nConvergence order with alpha =", alpha_val)
    print(calculate_order_of_convergence(actual_derivative, forward_diff_results, alpha_val))
    print("\n")

    # Centered difference method
    print("CENTERED DIFFERENCE METHOD")
    centered_diff_results = calculate_centered_difference(func, h_values, x_val)
    print("Approximations:")
    print(centered_diff_results)
    alpha_val = 1
    print("\nConvergence order with alpha =", alpha_val)
    print(calculate_order_of_convergence(actual_derivative, centered_diff_results, alpha_val))
    print("\n")

    print("Difference between methods / h:")
    print(np.abs(forward_diff_results - centered_diff_results))
    return


def calculate_forward_difference(func, h_values, x_val):
    func_at_x = func(x_val)
    approximations = np.zeros(len(h_values))

    for i, h in enumerate(h_values):
        approximations[i] = func(x_val + h) - func_at_x

    return approximations / h_values


def calculate_centered_difference(func, h_values, x_val):
    approximations = np.zeros(len(h_values))

    for i, h in enumerate(h_values):
        approximations[i] = func(x_val + h) - func(x_val - h)

    return approximations / (2 * h_values)


def calculate_order_of_convergence(actual, approximations, alpha):
    seq = np.zeros(len(approximations) - 1)

    for n, approx in enumerate(approximations):
        if n == len(approximations) - 1:
            continue

        seq[n] = np.abs(approximations[n + 1] - actual) / (np.abs(approx - actual)) ** alpha

    return seq


print("\n")
calculate_derivatives()
print("\n")
