import matplotlib
import matplotlib.pyplot as plt
import math
from math import cos, sin
import numpy as np
from matplotlib.ticker import LinearLocator
from matplotlib import cm
matplotlib.use('TkAgg')


def theta_dot(u):
    return u

def phi_dot(v):
    return v

def psi_dot(w):
    return w

def u_dot(theta, v, w):
    v_2 = v**2
    sin_cos = sin(theta) * cos(theta)
    return (
        (v_2 * sin_cos) - ((C/A) * ((w * v * sin(theta))) + (v_2 * sin_cos)) +
        ((M * 9.81 * h / A) * sin(theta))
    )

def v_dot(theta, u, v, w):
    sin_cos = sin(theta) * cos(theta)
    return (
        (
            ((C/A) * ((w * u * sin(theta)) + (v * u * sin_cos)))
            - (2 * v * u * sin_cos)
        )
        / sin(theta)**2
    )

def w_dot(theta, u, v, w):
    sin_cos = sin(theta) * cos(theta)
    return (v * u * sin(theta)) - (
        (
            (
                ((C/A) * ((w * u * sin(theta)) + (v * u * sin_cos)))
                - (2 * v * u * sin_cos)
            ) * cos(theta)
        )
        / sin(theta)**2
    )

def runge_kutta(theta1, phi1, psi1, u1, v1, w1):
    h = 0.001

    steps = int(abs((2 * math.pi) - phi1) // h)
    
    theta = [theta1]
    phi = [phi1]
    psi = [psi1]
    u = [u1]
    v = [v1]
    w = [w1]
    for _ in range(steps):
        c1 = theta_dot(u[-1])
        d1 = phi_dot(v[-1])
        e1 = psi_dot(w[-1])
        f1 = u_dot(theta[-1], v[-1], w[-1])
        g1 = v_dot(theta[-1], u[-1], v[-1], w[-1])
        i1 = w_dot(theta[-1], u[-1], v[-1], w[-1])

        c2 = theta_dot(u[-1] + (h * f1 / 2))
        d2 = phi_dot(v[-1] + (h * g1 / 2))
        e2 = psi_dot(w[-1] + (h * i1 / 2))
        f2 = u_dot(
            theta[-1] + (h * c1 / 2),
            v[-1] + (h * g1 / 2),
            w[-1] + (h * i1 / 2)
        )
        g2 = v_dot(
            theta[-1] + (h * c1 / 2),
            u[-1] + (h * f1 / 2),
            v[-1] + (h * g1 / 2),
            w[-1] + (h * i1 / 2)
        )
        i2 = w_dot(
            theta[-1] + (h * c1 / 2),
            u[-1] + (h * f1 / 2),
            v[-1] + (h * g1 / 2),
            w[-1] + (h * i1 / 2)
        )

        c3 = theta_dot(u[-1] + (h * f2 / 2))
        d3 = phi_dot(v[-1] + (h * g2 / 2))
        e3 = psi_dot(w[-1] + (h * i2 / 2))
        f3 = u_dot(
            theta[-1] + (h * c2 / 2),
            v[-1] + (h * g2 / 2),
            w[-1] + (h * i2 / 2)
        )
        g3 = v_dot(
            theta[-1] + (h * c2 / 2),
            u[-1] + (h * f2 / 2),
            v[-1] + (h * g2 / 2),
            w[-1] + (h * i2 / 2)
        )
        i3 = w_dot(
            theta[-1] + (h * c2 / 2),
            u[-1] + (h * f2 / 2),
            v[-1] + (h * g2 / 2),
            w[-1] + (h * i2 / 2)
        )

        c4 = theta_dot(u[-1] + (h * f3))
        d4 = phi_dot(v[-1] + (h * g3))
        e4 = psi_dot(w[-1] + (h * i3))
        f4 = u_dot(
            theta[-1] + (h * c3),
            v[-1] + (h * g3),
            w[-1] + (h * i3)
        )
        g4 = v_dot(
            theta[-1] + (h * c3),
            u[-1] + (h * f3),
            v[-1] + (h * g3),
            w[-1] + (h * i3)
        )
        i4 = w_dot(
            theta[-1] + (h * c3),
            u[-1] + (h * f3),
            v[-1] + (h * g3),
            w[-1] + (h * i3)
        )

        theta.append(theta[-1] + ((h / 6) * (c1 + (2 * c2) + (2 * c3) + c4)))
        phi.append(phi[-1] + ((h / 6) * (d1 + (2 * d2) + (2 * d3) + d4)))
        psi.append(psi[-1] + ((h / 6) * (e1 + (2 * e2) + (2 * e3) + e4)))
        u.append(u[-1] + ((h / 6) * (f1 + (2 * f2) + (2 * f3) + f4)))
        v.append(v[-1] + ((h / 6) * (g1 + (2 * g2) + (2 * g3) + g4)))
        w.append(w[-1] + ((h / 6) * (i1 + (2 * i2) + (2 * i3) + i4)))

    return np.array(theta), np.array(phi), np.array(psi)

C = 32/10000
A = 25/10000
M = 1
h = 3/100

def main():
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})

    theta1 = math.pi / 3
    phi1 = 0
    psi1 = 0
    u1 = 0
    omega = 78 # 200, 78, 5, 1.4, 0.7, 0
    w1 = 20 * math.pi
    
    theta, phi, psi = runge_kutta(theta1, phi1, psi1, u1, omega, w1)
    x = np.sin(theta) * np.sin(phi)
    y = np.sin(theta) * np.cos(phi)
    z = np.cos(theta)

    # Plot the trajectory of the symmetry axis
    ax.plot(x, y, z, color="black")

    # Plot the sphere for visual reference
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(
        sphere_x, sphere_y, sphere_z, color='b', alpha=0.1, edgecolor='k'
    )

    # Set plot labels and limits
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(
        f"Trajectory of the symmetry axis on a sphere for Omega: {omega}"
    )
    ax.legend()

    plt.show()

if __name__ == "__main__":
    main()
