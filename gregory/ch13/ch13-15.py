import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.ticker import LinearLocator
from matplotlib import cm
matplotlib.use('TkAgg')


def dy_fn(u):
    return u

def du_fn(x0, x, u):
    u_2 = u**2
    four_x0_2 = (4*(x0**2))
    eight_x_2 = (8*(x**2))
    x0_2 = x0**2
    x_2 = x**2
    four_x_2 = (4*(x**2))

    return (
        x * u * (four_x0_2 - eight_x_2 - 1 - u_2)
    ) / ((x0_2 - x_2) * (1 + four_x_2))

def runge_kutta(x0, x1, y1, u1):
    h = 0.01

    steps = int(abs(x1 - x0) // h)
    
    y = [y1]
    u = [u1]
    x = [x1]
    for _ in range(steps):
        c1 = dy_fn(u[-1])
        d1 = du_fn(x0, x[-1], u[-1])

        c2 = dy_fn(u[-1] + (h * d1 / 2))
        d2 = du_fn(x0, x[-1] + (h / 2), u[-1] + (h * d1 / 2))

        c3 = dy_fn(u[-1] + (h * d2 / 2))
        d3 = du_fn(x0, x[-1] + (h / 2), u[-1] + (h * d2 / 2))

        c4 = dy_fn(u[-1] + (h * d3))
        d4 = du_fn(x0, x[-1] + h, u[-1] + (h * d3))

        y.append(y[-1] + ((h / 6) * (c1 + (2 * c2) + (2 * c3) + c4)))
        u.append(u[-1] + ((h / 6) * (d1 + (2 * d2) + (2 * d3) + d4)))
        x.append(x[-1] + h)

    return np.array(x), np.array(y) 


def main():
    fig, ax1 = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})
    fig, ax2 = plt.subplots(figsize=(10, 10))

    x0 = 2
    x1 = 0
    y1 = 1

    legends = []
    for _lambda in [-0.2, -0.3, -0.31, -0.6]:
        x, y = runge_kutta(x0, x1, y1, u1=_lambda)
        ax2.plot(x, y)
        legends.append(f"lambda: {_lambda}")

    ax2.legend(labels=legends)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid(True)

    # 3D plot
    x, y = runge_kutta(x0, x1, y1, u1=-0.31)
    X = np.arange(-3, 3, 0.25)
    Y = np.arange(-2, 2, 0.25)
    X, Y = np.meshgrid(X, Y)
    Z = X**2
    ax1.plot(x, y, color="red", linestyle="dashed")
    ax1.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
    # surf = ax1.plot_surface(X, Y, Z, linewidth=0, antialiased=False)
    ax1.scatter(2, 0, 4, c='black', marker='o', s=100)
    ax1.scatter(0, 1, 0, c='black', marker='o', s=100)

    ax1.scatter(x, y, x**2, c='red', marker='.')

    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')

    plt.show()


if __name__ == "__main__":
    main()
