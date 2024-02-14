import numpy as np
from mod_fixedpt import find_fixed_point

def aitkens_acceleration(sequence, sequence_length, tolerance):
    accelerated_sequence = []
    for (value, i) in list(zip(sequence, range(sequence_length))):
        new_value = sequence[i] - (((sequence[i + 1] - sequence[i]) ** 2) / (sequence[i + 2] - (2 * sequence[i + 1]) + sequence[i]))
        accelerated_sequence.append(new_value)
        if (abs(accelerated_sequence[i] - accelerated_sequence[i - 1]) / abs(accelerated_sequence[i - 1]) < tolerance and i != 0):
            return (accelerated_sequence, i)
    
    return accelerated_sequence

# Define the function for which fixed point is to be found
function = lambda x: ((10) / (x + 4)) ** (1 / 2)

# Apply fixed point iteration method
approximate_fixed_point, error_code, iterations = find_fixed_point(function, initial_guess=1.5, tolerance=1e-10, max_iterations=1000)
print("Fixed point approximation: " + str(approximate_fixed_point[-1]))  # Use -1 to get the last element
print("Iterations needed: " + str(iterations))
print("Error code: " + str(error_code))

# Apply Aitken's method for acceleration
accelerated_fixed_point, count = aitkens_acceleration(approximate_fixed_point, iterations, tolerance=1e-10)
print("")
print("Aitken's Method Approximation: " + str(accelerated_fixed_point[count]))
print("Iterations needed: " + str(count))
