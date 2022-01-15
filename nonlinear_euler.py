import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    plt.figure()
    ax = plt.axes()
    t, x = nonlinear_euler(x0=1, t0=0, h=0.02, tm=3000)
    plot(ax, t, x, "Euler approximation h=0.02 and tm=3000", )
    t, x = nonlinear_euler(x0=1, t0=0, h=0.02, tm=900)
    plot(ax, t, x, "Euler approximation h=0.02 and tm=900", )
    plt.savefig("nonlinear_euler.png")
    plt.show()


def nonlinear_euler(x0, t0, h=0.05, tm=1):
    x, t = [x0], [t0]

    while abs(t[-1] - tm) > (h / 2):
        x.append(x[-1] + (h * (t[-1] - x[-1]**2)))
        t.append(t[-1] + h)
    return t, x


def plot(ax, t, x, label):
    ax.plot(t, x, label=label)
    ax.legend()
    ax.set_xlabel("t")
    ax.set_ylabel("x")


if __name__ == "__main__":
    main()
