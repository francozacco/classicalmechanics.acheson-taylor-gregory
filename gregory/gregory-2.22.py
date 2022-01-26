import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    plt.figure()
    ax = plt.axes()
    x, y = runge_kutta(x0=1, y0=1, vd=3, vhx=3, vhy=0)
    plot(ax, x, y)
    plt.show()


# def imprvd_linear_euler(x0, t0, h=0.05, tm=1):
#     x, t = x0, t0
#     while abs(t - tm) > (h / 2):
#         c1 = h * x
#         c2 = h * (x + c1)
#         x = x + 0.5 * (c1 + c2)
#         t = t + h
#     print(f"h={h}, x(1)={x}, error={math.e - x}")


def X(x, y, vd, vh):
    return ((-vd*x) / math.sqrt(x**2 + y**2)) - vh


def Y(x, y, vd, vh):
    return ((-vd*y) / math.sqrt(x**2 + y**2)) - vh


def runge_kutta(x0, y0, vd, vhx, vhy, h=0.05, xm=10):
    x, y = [x0], [y0]
    for _ in range(xm):
        c1 = h * X(x[-1], y[-1], vd, vhx)
        d1 = h * Y(x[-1], y[-1], vd, vhy)
        c2 = h * X(x[-1] + (h/2)*c1, y[-1] + (h/2)*d1, vd, vh)
        d2 = h * Y(x[-1] + (h/2)*c1, y[-1] + (h/2)*d1, vd, vh)
        c3 = h * X(x[-1] + (h/2)*c2, y[-1] + (h/2)*d2, vd, vh)
        d3 = h * Y(x[-1] + (h/2)*c2, y[-1] + (h/2)*d2, vd, vh)
        c4 = h * X(x[-1] + c3, y[-1] + d3, vd, vh)
        d4 = h * Y(x[-1] + c3, y[-1] + d3, vd, vh)
        x.append(x[-1] + 1/6 * (c1 + (2 * c2) + (2 * c3) + c4))
        y.append(y[-1] + 1/6 * (d1 + (2 * d2) + (2 * d3) + d4))
    return x, y
    # print(f"h={h}, x(1)={x}, error={math.e - x}")


def plot(ax, x, y):
    print(x,y)
    ax.plot(x, y)  # label=)
    # ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")


if __name__ == "__main__":
    main()
