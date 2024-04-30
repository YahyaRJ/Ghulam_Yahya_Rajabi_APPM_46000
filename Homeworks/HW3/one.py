import numpy as np

# import bisection method from class
# with added print statements at each iteration
# and tolerance checking for relative error
# instead of absolute error


def driver():
    # function decleration
    f = lambda x: 2 * x - 1 - np.sin(x)
    

    # endpoints for parts a,b,c
    a = 0
    b = np.pi / 2

    # tolerance for 8 correct digits
    tol = 0.5 * 10**-8

    print("(1):\n")
    [astar, ier] = bisection(f, a, b, tol)
    print("the approximate root is", astar)
    # print("the error message reads:", ier)
    print("f(root) =", f(astar))
    print("\n")

    return


# define routines
def bisection(f, a, b, tol):
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        ier = 1
        astar = a
        return [astar, ier]

    #   verify end points are not a root
    if fa == 0:
        astar = a
        ier = 0
        return [astar, ier]

    if fb == 0:
        astar = b
        ier = 0
        return [astar, ier]

    count = 0
    d = 0.5 * (a + b)
    while abs(d - a) / abs(d) > tol:
        fd = f(d)
        if fd == 0:
            astar = d
            ier = 0
            return [astar, ier]
        if fa * fd < 0:
            b = d
        else:
            a = d
            fa = fd
        d = 0.5 * (a + b)
        count = count + 1
        print("iteration: ", count, " | curr_root = ", d)
    #      print('abs(d-a) = ', abs(d-a))

    print("Number of iterations: ", count, "\n")
    astar = d
    ier = 0
    return [astar, ier]


print("")
driver()