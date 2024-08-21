import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np


def main():
    plt.figure()
    ax = plt.axes()

    xi_min, xi_max, xi_step = 0, 1, 0.02
    F = []
    xi_range = np.arange(xi_min, xi_max, xi_step)
    for xi in xi_range:
        F.append(
            integrate(xi=xi, theta_min=0, theta_max=2*math.pi)
        )
    print(xi_range[0], F[0])
    plot(ax, xi_range, F)
    plt.show()
    # plt.savefig("gregory-3.9.png")


def f(theta, xi):
    return (1 - xi*math.cos(theta)) / (1 + xi**2 - 2*xi*math.cos(theta))**(3/2)


def integrate(xi, theta_min, theta_max, steps=1000):
    h = (theta_max - theta_min) / steps
    x = theta_min
    res = f(x, xi)
    res += f(theta_max, xi)
    for i in range(1, steps):
        x += h
        res += 2 * f(x, xi) if i % 2 == 0 else 4 * f(x, xi)

    return (h / 3) * res


def plot(ax, x, y):
    ax.plot(x, y)
    ax.set_xlabel("\u03BE")
    ax.set_ylabel("F")


if __name__ == "__main__":
    main()
