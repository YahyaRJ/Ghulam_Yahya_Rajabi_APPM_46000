import numpy as np

def main():
    # Define test functions
    # The true fixed point for f1 is approximately 1.4987
    function1 = lambda x: 1 + 0.5 * np.sin(x)
    # The true fixed point for f2 is approximately 3.09
    function2 = lambda x: 3 + 2 * np.sin(x)

    # Set maximum iterations and tolerance
    max_iterations = 100
    tolerance = 1e-6

    # Test function 1
    initial_guess = 0.0
    approximate_fixed_point, error_code, iterations = find_fixed_point(function1, initial_guess, tolerance, max_iterations)
    print('Approximate fixed point:', approximate_fixed_point[iterations])
    print('Value of f1(approximate fixed point):', function1(approximate_fixed_point[iterations]))
    print('Error code:', error_code)

    # Test function 2
    initial_guess = 0.0
    approximate_fixed_point, error_code, iterations = find_fixed_point(function2, initial_guess, tolerance, max_iterations)
    print('Approximate fixed point:', approximate_fixed_point[iterations])
    print('Value of f2(approximate fixed point):', function2(approximate_fixed_point[iterations]))
    print('Error code:', error_code)

def find_fixed_point(f, initial_guess, tolerance, max_iterations):
    ''' 
    Finds the fixed point of a function using the fixed-point iteration method.

    Parameters: 
    f (function): The function for which the fixed point is to be found.
    initial_guess (float): The initial guess for the fixed point.
    tolerance (float): The tolerance for stopping iterations.
    max_iterations (int): The maximum number of iterations allowed.

    Returns: 
    approximate_fixed_point (float): The approximate fixed point.
    error_code (int): 0 for success, 1 for error.
    iterations (int): The number of iterations performed.
    '''

    # Initialize iteration counter
    iteration_count = 0

    # Initialize array to store all approximations
    approximation_history = np.zeros(max_iterations + 1)
    # Set initial value
    approximation_history[iteration_count] = initial_guess

    # Iterate until maximum iterations reached
    while iteration_count < max_iterations:
        # Increment iteration counter
        iteration_count += 1

        # Calculate next approximation using fixed-point iteration formula: x_(n+1) = f(x_n)
        next_approximation = f(initial_guess)
        # Store the new approximation
        approximation_history[iteration_count] = next_approximation

        # Check if tolerance criteria met
        if abs(next_approximation - initial_guess) / abs(initial_guess) < tolerance:
            approximate_fixed_point = next_approximation
            # Success
            error_code = 0
            return approximation_history, error_code, iteration_count

        # Update initial guess for the next iteration
        initial_guess = next_approximation

    # If maximum iterations reached without meeting tolerance criteria
    approximate_fixed_point = next_approximation
    # Failure
    error_code = 1
    return approximation_history, error_code, iteration_count

#main()
