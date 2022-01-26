import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    plt.figure()
    ax = plt.axes()
    t, x = linear_euler(x0=1, t0=0)
    plot(ax, t, x, "Euler approximation")
    plot(ax, t, [math.e**i for i in t], "Exact solution")
    plt.savefig("euler.png")
    plt.show()


def linear_euler(x0, t0, h=0.05, tm=1):
    x, t = [x0], [t0]
    while abs(t[-1] - tm) > (h / 2):
        x.append(x[-1] + (h * x[-1]))
        t.append(t[-1] + h)
    print("Result: ", x[-1])
    return t, x


def plot(ax, t, x, label):
    ax.plot(t, x, label=label)
    ax.legend()
    ax.set_xlabel("t")
    ax.set_ylabel("x")


if __name__ == "__main__":
    main()
