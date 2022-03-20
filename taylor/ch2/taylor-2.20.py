import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np


def main():
    plt.figure()
    ax = plt.axes()

    t = np.arange(0, 3, 0.01)
    for i in [0.3, 1, 3]:
        calculate_curve(ax, v_ter=i, t=t)

    # v_ter = inf case:
    plot(ax, t, [ti - ti**2/2 for ti in t], "v_ter = inf")
    plt.show()
    # plt.savefig("taylor-2.20.png")



def calculate_curve(ax, v_ter, t):
    x = []
    y = []
    tau = v_ter

    for ti in t:
        x_f = tau * (1 - math.e**(-ti / tau))
        x.append(x_f)
        y_f = (1 + v_ter) * tau * (1 - math.e**(-ti / tau)) - v_ter * ti
        y.append(y_f)

    plot(ax, x, y, f"v_ter = {v_ter}")


def plot(ax, x, y, label):
    ax.plot(x, y, label=label)
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")


if __name__ == "__main__":
    main()
