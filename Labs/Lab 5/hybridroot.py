# Import libraries
import numpy as np

# Define the main driver function
def driver():
    # Define the functions f, f', and f'' for the equation
    function = lambda x: np.exp((x ** 2) + (7 * x) - 30) - 1
    derivative = lambda x: ((2*x) + 7) * np.exp((x ** 2) + (7 * x) - 30)
    second_derivative = lambda x: (((2*x) + 7) ** 2) * np.exp((x ** 2) + (7 * x) - 30)
    
    # Define the initial interval [left_endpoint, right_endpoint], tolerance, and maximum number of iterations
    left_endpoint = 2
    right_endpoint = 4.5
    tolerance = 1e-14
    max_iterations = 100

    # Call the hybrid function to find the root using the hybrid method
    [approx_root, error_code, iterations] = hybrid(function, derivative, second_derivative, left_endpoint, right_endpoint, max_iterations, tolerance)

    # Print the results
    print('The number of iterations was ', iterations)
    print('The approximate root is', approx_root)
    print('The error message reads:', error_code)
    print('f(pstar) =', function(approx_root))

# Define the hybrid root-finding function
def hybrid(function, derivative, second_derivative, left_endpoint, right_endpoint, max_iterations, tolerance):
    # Evaluate f(left_endpoint) and f(right_endpoint)
    f_left = function(left_endpoint)
    f_right = function(right_endpoint)
    
    # Check if the signs of f(left_endpoint) and f(right_endpoint) are the same, indicating no root in the interval
    if (f_left * f_right > 0):
        error_code = 1
        root_approximation = left_endpoint
        return [root_approximation, error_code]

    # Check if the endpoints are roots
    if (f_left == 0):
        root_approximation = left_endpoint
        error_code = 0
        return [root_approximation, error_code, 0]

    if (f_right == 0):
        root_approximation = right_endpoint
        error_code = 0
        return [root_approximation, error_code, 0]

    # Initialize iteration count
    iterations = 0

    # Initial guess using midpoint of the interval
    guess = 0.5 * (left_endpoint + right_endpoint)
    
    # Compute the ratio gp which is used to check convergence
    convergence_ratio = abs((function(guess) * second_derivative(guess)) / (derivative(guess) ** 2))
    
    # Iteratively refine the interval until convergence
    while (convergence_ratio > 1):
        f_guess = function(guess)
        if (f_guess == 0):
            root_approximation = guess
            error_code = 0
            return [root_approximation, error_code, iterations]
        if (f_left * f_guess < 0):
            right_endpoint = guess
        else:
            left_endpoint = guess
            f_left = f_guess
        guess = 0.5 * (left_endpoint + right_endpoint)
        convergence_ratio = abs((function(guess) * second_derivative(guess)) / (derivative(guess) ** 2))
        iterations += 1

    # Use the initial guess as the first iteration for Newton's method
    initial_approximation = guess
    approximations = np.zeros(max_iterations + 1)
    approximations[0] = guess
    
    # Iterate using Newton's method
    for iteration in range(max_iterations):
        iterations += 1
        next_approximation = initial_approximation - function(initial_approximation) / derivative(initial_approximation)
        approximations[iteration + 1] = next_approximation
        if (abs(next_approximation - initial_approximation) < tolerance):
            final_root = next_approximation
            error_code = 0
            return [final_root, error_code, iterations]
        initial_approximation = next_approximation

    # If maximum iterations reached without convergence, return the last approximation
    final_root = next_approximation
    error_code = 1
    return [final_root, error_code, iterations]

# Call the driver function to execute the program
driver()
