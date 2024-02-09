import numpy as np
import matplotlib.pyplot as plt
import random
import math

# Calculate condition number
matrix_A = 0.5 * np.array([[1, 1], [1 + 10 ** -10, 1 - 10 ** -10]])
poi_matrix = np.array([[1 - 10 ** 10, 10 ** 10], [1 + 10 ** 10, -1 * 10 ** 10]])
max_singular_val = np.linalg.svd(matrix_A)[1][0]
min_singular_val = np.linalg.svd(matrix_A)[1][1]
condition_num = max_singular_val / min_singular_val
print("Condition number: " + str(condition_num))

# Calculate norm and relative error
delta_B = np.array([5 * 10 ** -5, 5 * 10 ** -5])
delta_X = np.matmul(poi_matrix, delta_B)
norm_X = np.linalg.norm(delta_X)
print("Norm: " + str(norm_X))
relative_error = condition_num * norm_X
print("Relative error: " + str(relative_error))

# Define a function
def exponential_function(x):
    y = (math.e) ** x
    return y - 1

# Test the function
input_value = np.longdouble(9.999999995000000 * 10 ** -10)
print("Function output: " + str(exponential_function(input_value)))

# Calculate relative error iteratively
def calculate_relative_error():
    # Initialize variables
    relative_error_val = np.longdouble(10 ** 10)
    actual_value = np.longdouble(10 ** -9)
    accuracy_tolerance = np.longdouble(10 ** -16)
    g_n = np.longdouble(0)
    j = 0
    while (relative_error_val > accuracy_tolerance):
        j += 1
        g_n += np.longdouble((input_value ** j)) / np.longdouble((math.factorial(j)))
        relative_error_val = np.absolute(actual_value - g_n) / actual_value

    print("Number of terms: " + str(j))
    return relative_error_val

relative_error_input_value = calculate_relative_error()

# Define a function
def g_2_function(x):
    return x + ((x ** 2) / 2)

# Print the function output and relative error
print("g_2 function output: " + "{0:.40f}".format(g_2_function(np.longdouble(input_value))))
print("Relative error: " + "{0:.40f}".format(relative_error_input_value))

# Calculate the sum
start_val = 0
end_val = np.pi
interval_points = np.linspace(start_val, end_val, 31) 
cos_values = np.cos(interval_points)
sum_result = 0
for point, cos_val in zip(interval_points, cos_values):
    sum_result += point * cos_val
print("Sum result: " + str(sum_result))

# Plot curves
theta_values = np.linspace(0, 2 * np.pi, 1000)
amplitude = 1.2
delta_radius = 0.1
frequency = 15
phase = 0
plt.figure()
plt.plot()

plt.figure()
delta_radius = 0.05
for i in range(1, 11):
    amplitude_val = i
    frequency_val = 2 + i
    phase_val = random.uniform(0, 2)
    x_curve = amplitude_val * (1 + delta_radius * np.sin(frequency_val * theta_values + phase_val)) * np.cos(theta_values)
    y_curve = amplitude_val * (1 + delta_radius * np.sin(frequency_val * theta_values + phase_val)) * np.sin(theta_values)
    plt.plot(x_curve, y_curve)
plt.show()
