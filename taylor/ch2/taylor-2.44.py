import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math

D = 0.15
# D = 0.
dens = 7800
v0 = 300
m = 13.78
gamma = 0.25
lambd = 10000
theta = 50 * math.pi / 180


def vacuum_traj(x0, y0, v0, h, it):
    x, y = [x0], [y0]
    vx0, vy0 = v0*math.cos(theta), v0*math.sin(theta)
    t = 0

    for _ in range(it):
        x.append((vx0*t))
        y.append((vy0*t) - (4.9*(t**2)))
        t = t + h

    return x, y


def fc(y):
    return gamma * (D**2) * (math.e**(-y/lambd))


def fx(vx):
    return vx


def fy(vy):
    return vy


def fvx(x, y, vx, vy, c):
    return - (c(y)/m) * math.sqrt(vx**2 + vy**2) * vx


def fvy(x, y, vx, vy, c):
    return - 9.8 - (c(y)/m) * math.sqrt(vx**2 + vy**2) * vy


def runge_kutta(x0, y0, v0, c=fc, h=0.01, it=1500):
    x, y = [x0], [y0]
    vx, vy = [v0*math.cos(theta)], [v0*math.sin(theta)]
    t = 0

    for _ in range(it):
        c1 = h * fvx(x[-1], y[-1], vx[-1], vy[-1], c)
        d1 = h * fvy(x[-1], y[-1], vx[-1], vy[-1], c)
        e1 = h * fx(vx[-1])
        f1 = h * fy(vy[-1])

        c2 = h * fvx(
            x[-1] + (h/2)*e1, y[-1] + (h/2)*f1,
            vx[-1] + (h/2)*c1, vy[-1] + (h/2)*d1,
            c
        )
        d2 = h * fvy(
            x[-1] + (h/2)*e1, y[-1] + (h/2)*f1,
            vx[-1] + (h/2)*c1, vy[-1] + (h/2)*d1,
            c
        )
        e2 = h * fx(vx[-1] + (h/2)*c1)
        f2 = h * fy(vy[-1] + (h/2)*d1)

        c3 = h * fvx(
            x[-1] + (h/2)*e2, y[-1] + (h/2)*f2,
            vx[-1] + (h/2)*c2, vy[-1] + (h/2)*d2,
            c
        )
        d3 = h * fvy(
            x[-1] + (h/2)*e2, y[-1] + (h/2)*f2,
            vx[-1] + (h/2)*c2, vy[-1] + (h/2)*d2,
            c
        )
        e3 = h * fx(vx[-1] + (h/2)*c2)
        f3 = h * fy(vy[-1] + (h/2)*d2)

        c4 = h * fvx(x[-1] + e3, y[-1] + f3, vx[-1] + c3, vy[-1] + d3, c)
        d4 = h * fvy(x[-1] + e3, y[-1] + f3, vx[-1] + c3, vy[-1] + d3, c)
        e4 = h * fx(vx[-1] + c3)
        f4 = h * fy(vy[-1] + d3)

        vx.append(vx[-1] + 1/6 * (c1 + (2 * c2) + (2 * c3) + c4))
        vy.append(vy[-1] + 1/6 * (d1 + (2 * d2) + (2 * d3) + d4))
        x.append(x[-1] + 1/6 * (e1 + (2 * e2) + (2 * e3) + e4))
        y.append(y[-1] + 1/6 * (f1 + (2 * f2) + (2 * f3) + f4))
        t = t + h

    return x, y, vx, vy


def plot(ax, x, y, label):
    ax.plot(x, y, label=label)
    ax.legend()
    ax.set_ylim(bottom=0, top=2750)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def main():
    plt.figure()
    ax = plt.axes()
    t = 35
    h = 0.01
    it = int(t/h)

    x, y, vx, vy = runge_kutta(x0=0, y0=0, v0=v0, c=fc, h=h, it=it)
    plot(ax, x, y, "Trajectory - c(y)")
    x, y, vx, vy = runge_kutta(x0=0, y0=0, v0=v0, c=lambda x: fc(0), h=h, it=it)
    plot(ax, x, y, "Trajectory - c(0)")
    x, y = vacuum_traj(x0=0, y0=0, v0=v0, h=h, it=it)
    plot(ax, x, y, "Trajectory - vacuum")

    plt.show()
    # plt.savefig("taylor-2.44.png")


if __name__ == "__main__":
    main()
