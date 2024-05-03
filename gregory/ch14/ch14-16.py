import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np

matplotlib.use('TkAgg')


def x_fn(p_x, m):
    return p_x / m

def y_fn(p_y, m):
    return p_y / m

def p_x_fn(x, y, t, m):
    x_m_cos = x - (math.cos(t))
    x_m_sq = (x_m_cos)**2
    y_m_sin = y - (math.sin(t))
    y_m_sq = (y_m_sin)**2
    x_p_cos = x + (math.cos(t))
    x_p_sq = (x_p_cos)**2
    y_p_sin = y + (math.sin(t))
    y_p_sq = (y_p_sin)**2
    return -4 * m * (
        (x_m_cos / ((x_m_sq + y_m_sq)**(3/2)))
        + (x_p_cos / ((x_p_sq + y_p_sq)**(3/2)))
    )

def p_y_fn(x, y, t, m):
    x_m_cos = x - (math.cos(t))
    x_m_sq = (x_m_cos)**2
    y_m_sin = y - (math.sin(t))
    y_m_sq = (y_m_sin)**2
    x_p_cos = x + (math.cos(t))
    x_p_sq = (x_p_cos)**2
    y_p_sin = y + (math.sin(t))
    y_p_sq = (y_p_sin)**2
    return -4 * m * (
        (y_m_sin / ((x_m_sq + y_m_sq)**(3/2)))
        + (y_p_sin / ((x_p_sq + y_p_sq)**(3/2)))
    )


def runge_kutta(
    ax1, m=1, x_init=1, y_init=1, p_x_init=1, p_y_init=1,
):
    h = 0.001
    steps = 100000

    x = [x_init]
    y = [y_init]
    p_x = [p_x_init]
    p_y = [p_y_init]
    t = [0]
    for _ in range(steps):
        c1 = x_fn(p_x[-1], m)
        d1 = y_fn(p_y[-1], m)
        e1 = p_x_fn(x[-1], y[-1], t[-1], m)
        f1 = p_y_fn(x[-1], y[-1], t[-1], m)

        c2 = x_fn(p_x[-1] + (h * e1 / 2), m)
        d2 = y_fn(p_y[-1] + (h * f1 / 2), m)
        e2 = p_x_fn(
            x[-1] + (h * c1 / 2), y[-1] + (h * d1 / 2), t[-1] + (h / 2), m
        )
        f2 = p_y_fn(
            x[-1] + (h * c1 / 2), y[-1] + (h * d1 / 2), t[-1] + (h / 2), m
        )

        c3 = x_fn(p_x[-1] + (h * e2 / 2), m)
        d3 = y_fn(p_y[-1] + (h * f2 / 2), m)
        e3 = p_x_fn(
            x[-1] + (h * c2 / 2), y[-1] + (h * d2 / 2), t[-1] + (h / 2), m
        )
        f3 = p_y_fn(
            x[-1] + (h * c2 / 2), y[-1] + (h * d2 / 2), t[-1] + (h / 2), m
        )

        c4 = x_fn(p_x[-1] + (h * e3), m)
        d4 = y_fn(p_y[-1] + (h * f3), m)
        e4 = p_x_fn(x[-1] + (h * c3), y[-1] + (h * d3), t[-1] + h, m)
        f4 = p_y_fn(x[-1] + (h * c3), y[-1] + (h * d3), t[-1] + h, m)

        x.append(x[-1] + (h / 6 * (c1 + (2 * c2) + (2 * c3) + c4)))
        y.append(y[-1] + (h / 6 * (d1 + (2 * d2) + (2 * d3) + d4)))
        p_x.append(p_x[-1] + (h / 6 * (e1 + (2 * e2) + (2 * e3) + e4)))
        p_y.append(p_y[-1] + (h / 6 * (f1 + (2 * f2) + (2 * f3) + f4)))
        t.append(t[-1] + h)

    t = np.array(t)
    rot_x = np.array(x) * np.cos(-t) - np.array(y) * np.sin(-t)
    rot_y = np.array(x) * np.sin(-t) + np.array(y) * np.cos(-t)
    ax1.plot(rot_x, rot_y, color="black")

    # Primaries
    ax1.plot(1, 0, "o", color="green")
    ax1.plot(-1, 0, "o", color="green")


def main():
    _, ax1 = plt.subplots(figsize=(10, 10))
    legends = []
    m = 0.01
    x_init = 0
    y_init = 0.1
    p_x_init = 0
    p_y_init = 0

    runge_kutta(ax1, m, x_init, y_init, p_x_init, p_y_init)
    legends.append(
        f"m:{m}  x:{x_init}  y:{y_init}  px:{p_x_init}  py:{p_y_init}"
    )
    legends.append("Primaries")
    ax1.legend(labels=legends)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
