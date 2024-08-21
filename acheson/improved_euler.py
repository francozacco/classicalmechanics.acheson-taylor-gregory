import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    print("Improved Euler method")
    for h in [0.1, 0.01, 0.001]:
        imprvd_linear_euler(x0=1, t0=0, h=h)
    print()
    print("Runge-Kutta method")
    for h in [0.1, 0.01, 0.001]:
        runge_kutta(x0=1, t0=0, h=h)


def imprvd_linear_euler(x0, t0, h=0.05, tm=1):
    x, t = x0, t0
    while abs(t - tm) > (h / 2):
        c1 = h * x
        c2 = h * (x + c1)
        x = x + 0.5 * (c1 + c2)
        t = t + h
    print(f"h={h}, x(1)={x}, error={math.e - x}")


def runge_kutta(x0, t0, h=0.05, tm=1):
    x, t = x0, t0
    while abs(t - tm) > (h / 2):
        c1 = h * x
        c2 = h * (x + (0.5 * c1))
        c3 = h * (x + (0.5 * c2))
        c4 = h * (x + c3)
        x = x + 1/6 * (c1 + (2 * c2) + (2 * c3) + c4)
        t = t + h
    print(f"h={h}, x(1)={x}, error={math.e - x}")


def plot(ax, t, x, label):
    ax.plot(t, x, label=label)
    ax.legend()
    ax.set_xlabel("t")
    ax.set_ylabel("x")


if __name__ == "__main__":
    main()
