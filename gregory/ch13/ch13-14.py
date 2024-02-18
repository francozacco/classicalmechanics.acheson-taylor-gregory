import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.ticker import LinearLocator
from matplotlib import cm
matplotlib.use('TkAgg')


def dr_fn(u):
    return u

def du_fn(r, u):
    four_r2 = (4*(r**2))
    four_u2 = (4*(u**2))
    r_2 = r**2
    u_2 = u**2
    u_4 = u**4
    return (
        ((four_u2 + 1)*(r_2 + ((1 + four_r2)*u_2)))
        - (4*(four_r2 + 1)*u_4)
        - ((four_r2 - 1)*u_2)
    ) / (r*(four_r2 + 1))

def runge_kutta(r_init, u_init, color=None):
    h = 0.01
    t_init = 0
    t_fin = math.pi/2

    steps = int((t_fin - t_init) // h)

    r = [r_init]
    u = [u_init]
    t = [t_init]
    for _ in range(steps):
        c1 = dr_fn(u[-1])
        d1 = du_fn(r[-1], u[-1])

        c2 = dr_fn(u[-1] + (h * d1 / 2))
        d2 = du_fn(r[-1] + (h * c1 / 2), u[-1] + (h * d1 / 2))

        c3 = dr_fn(u[-1] + (h * d2 / 2))
        d3 = du_fn(r[-1] + (h * c2 / 2), u[-1] + (h * d2 / 2))

        c4 = dr_fn(u[-1] + (h * d3))
        d4 = du_fn(r[-1] + (h * c3), u[-1] + (h * d3))

        r.append(r[-1] + ((h / 6) * (c1 + (2 * c2) + (2 * c3) + c4)))
        u.append(u[-1] + ((h / 6) * (d1 + (2 * d2) + (2 * d3) + d4)))
        t.append(t[-1] + h)

    x = np.array(r) * np.cos(t)
    y = np.array(r) * np.sin(t)
    _x = np.concatenate([np.sort(x), x])
    _y = np.concatenate([np.sort(-y), y])
    return _x, _y 


def main():
    fig, ax1 = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})
    fig, ax2 = plt.subplots(figsize=(10, 10))
    u_init = 0

    legends = []
    for _lambda in [0.3, 0.6, 0.63, 0.9, 1]:
        x, y = runge_kutta(r_init=_lambda, u_init=u_init)
        ax2.plot(x, y)
        legends.append(f"lambda: {_lambda}")

    ax2.legend(labels=legends)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid(True)

    # 3D plot
    x, y = runge_kutta(r_init=0.63, u_init=0)
    X = np.arange(-2, 2, 0.25)
    Y = np.arange(-2, 2, 0.25)
    X, Y = np.meshgrid(X, Y)
    Z = X**2 + Y**2
    ax1.plot(x, y, color="red", linestyle="dashed")
    ax1.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
    # surf = ax1.plot_surface(X, Y, Z, linewidth=0, antialiased=False)
    ax1.scatter(0, 1, 1, c='black', marker='o', s=100)
    ax1.scatter(0, -1, 1, c='black', marker='o', s=100)

    ax1.scatter(x, y, x**2 + y**2, c='red', marker='.')

    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')

    plt.show()


if __name__ == "__main__":
    main()
