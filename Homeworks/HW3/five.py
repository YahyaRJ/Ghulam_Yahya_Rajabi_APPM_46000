import numpy as np
import matplotlib.pyplot as plt


def driver():
    # function definition
    f = lambda x: x - 4 * np.sin(2 * x) - 3

    # plot over x-axis that includes all roots
    X = np.arange(-np.pi, 2 * np.pi, 0.001)
    fig, ax1 = plt.subplots()
    ax1.plot(X, f(X))
    plt.grid(True)
    plt.show()

    return


driver()