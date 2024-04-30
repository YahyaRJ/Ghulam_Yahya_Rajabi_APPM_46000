import numpy as np
import matplotlib.pyplot as plt

# imported fixed-point method from provided example


def driver():
    # function definition
    f = lambda x: -np.sin(2 * x) + 5 / 4 * x - 3 / 4

    # tolerance for relative error to be accurate to 10 decimal digits
    tol = 0.5 * 10**-10
    Nmax = 100

    # defining inital points
    x1, x2, x3, x4, x5 = -0.9, -0.4, 1.7, 3, 4.5

    # running method
    print("\n(5):\n")
    print("looking for root at x=-0.898 with x0 =", x1)
    [xstar, ier] = fixedpt(f, x1, tol, Nmax)
    print("the approximate fixed point is:", xstar)
    print("f(fixed_point):", f(xstar))
    print("Error message reads:", ier)
    print("\n")

    print("looking for root at x=-0.544 with x0 =", x2)
    [xstar, ier] = fixedpt(f, x2, tol, Nmax)
    print("the approximate fixed point is:", xstar)
    print("f(fixed_point):", f(xstar))
    print("Error message reads:", ier)
    print("\n")

    print("looking for root at x=1.732 with x0 =", x3)
    [xstar, ier] = fixedpt(f, x3, tol, Nmax)
    print("the approximate fixed point is:", xstar)
    print("f(fixed_point):", f(xstar))
    print("Error message reads:", ier)
    print("\n")

    print("looking for root at x=3.162 with x0 =", x4)
    [xstar, ier] = fixedpt(f, x4, tol, Nmax)
    print("the approximate fixed point is:", xstar)
    print("f(fixed_point):", f(xstar))
    print("Error message reads:", ier)
    print("\n")

    print("looking for root at x=4.518 with x0 =", x5)
    [xstar, ier] = fixedpt(f, x5, tol, Nmax)
    print("the approximate fixed point is:", xstar)
    print("f(fixed_point):", f(xstar))
    print("Error message reads:", ier)
    print("\n")

    return


# define routines
def fixedpt(f, x0, tol, Nmax):
    """x0 = initial guess"""
    """ Nmax = max number of iterations"""
    """ tol = stopping tolerance"""

    count = 0
    while count < Nmax:
        count = count + 1
        x1 = f(x0)
        if abs(x1 - x0) / abs(x0) < tol:
            xstar = x1
            ier = 0
            return [xstar, ier]
        x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]


driver()