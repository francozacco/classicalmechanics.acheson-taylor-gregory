import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    plt.figure()
    ax = plt.axes()
    vd, vh = 0.99, 1
    ax.set_title(f'vd={vd}, vh={vh}')
    xD, yD, xH, yH = runge_kutta(x0=-1, y0=0, vd=vd, vh=vh)
    plot(ax, xH, yH, "Hare trajectory")
    plot(ax, xD, yD, "Dog trajectory")
    ax.plot(xH[-1], yH[-1], marker="o", color="blue")
    ax.plot(xD[-1], yD[-1], marker="o", color="orange")

    plt.show()
    # plt.savefig("gregory-2.22-2_vd=vh.png")


def fx(x, y, t, vd, vh):
    return ((-vd*x) / math.sqrt(x**2 + y**2)) + vh*math.sin(t)


def fy(x, y, t, vd, vh):
    return ((-vd*y) / math.sqrt(x**2 + y**2)) - vh*math.cos(t)


def runge_kutta(x0, y0, vd, vh, h=0.01, it=10000):
    x, y = [x0], [y0]
    xD, yD = [], []
    xH, yH = [], []
    t = 0

    for _ in range(it):
        c1 = h * fx(x[-1], y[-1], t, vd, vh)
        d1 = h * fy(x[-1], y[-1], t, vd, vh)
        c2 = h * fx(x[-1] + (h/2)*c1, y[-1] + (h/2)*d1, t, vd, vh)
        d2 = h * fy(x[-1] + (h/2)*c1, y[-1] + (h/2)*d1, t, vd, vh)
        c3 = h * fx(x[-1] + (h/2)*c2, y[-1] + (h/2)*d2, t, vd, vh)
        d3 = h * fy(x[-1] + (h/2)*c2, y[-1] + (h/2)*d2, t, vd, vh)
        c4 = h * fx(x[-1] + c3, y[-1] + d3, t, vd, vh)
        d4 = h * fy(x[-1] + c3, y[-1] + d3, t, vd, vh)
        x.append(x[-1] + 1/6 * (c1 + (2 * c2) + (2 * c3) + c4))
        y.append(y[-1] + 1/6 * (d1 + (2 * d2) + (2 * d3) + d4))

        xH.append(math.cos(t))
        yH.append(math.sin(t))

        xD.append(xH[-1] + x[-1])
        yD.append(yH[-1] + y[-1])
        t = t + h

    return xD, yD, xH, yH


def plot(ax, x, y, label):
    ax.plot(x, y, label=label)
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")


if __name__ == "__main__":
    main()
