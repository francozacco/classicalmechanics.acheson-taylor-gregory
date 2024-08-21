import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')


def main():
    e = 0.8
    eta = 0.005
    v_init = 1 / (1 + e)
    w_init = 0
    h = 0.0001
    max_theta = 16 * math.pi
    steps = int(max_theta / h)
    theta_range = np.arange(0, max_theta, h)

    v = [v_init]
    r = [1 / v_init]
    w = [w_init]
    for _ in range(steps - 1):
        c1 = w_fn(w[-1])
        d1 = v_fn(v[-1], e, eta)

        c2 = w_fn(w[-1] + (h * c1 / 2))
        d2 = v_fn(v[-1] + (h * d1 / 2), e, eta)

        c3 = w_fn(w[-1] + (h * c2 / 2))
        d3 = v_fn(v[-1] + (h * d2 / 2), e, eta)

        c4 = w_fn(w[-1] + (c3 * h))
        d4 = v_fn(v[-1] + (d3 * h), e, eta)
        v.append(v[-1] + (h / 6 * (c1 + (2 * c2) + (2 * c3) + c4)))
        w.append(w[-1] + (h / 6 * (d1 + (2 * d2) + (2 * d3) + d4)))
        r.append(1 / v[-1])

    y = [r[i] * math.sin(theta_range[i]) for i in range(len(v))]
    x = [r[i] * math.cos(theta_range[i]) for i in range(len(v))]

    plot(x, y)
    plt.show()
    # plt.savefig("gregory-7.28.png")


def v_fn(v, e, eta):
    return (1 / (1 - e**2)) + (eta * v**2) - v


def w_fn(w):
    return w


def plot(x, y):
    plt.figure()
    ax = plt.axes()
    ax.plot(x, y)
    ax.scatter([0], [0], marker="*", s=200, color="orange")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Mercury Perihelion Precession")


if __name__ == "__main__":
    main()
