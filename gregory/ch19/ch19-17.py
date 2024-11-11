import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
matplotlib.use('TkAgg')


def Lx_dot(L_y, L_z):
    return ((B - C) / (B * C)) * L_y * L_z

def Ly_dot(L_x, L_z):
    return ((C - A) / (A * C)) * L_x * L_z

def Lz_dot(L_x, L_y):
    return ((A - B) / (A * B)) * L_x * L_y

###############################################################################

def e1x_dot():
    return 0

def e1y_dot(L_z):
    return L_z / C

def e1z_dot(L_y):
    return -(L_y / B)

###############################################################################

def e2x_dot(L_z):
    return -(L_z / C)

def e2y_dot():
    return 0

def e2z_dot(L_x):
    return L_x / A

###############################################################################

def e3x_dot(L_y):
    return L_y / B

def e3y_dot(L_x):
    return -(L_x / A)

def e3z_dot():
    return 0


def runge_kutta(L_x1, L_y1, L_z1):
    h = 0.01

    steps = int(100 // h)
    
    L_x, L_y, L_z = [L_x1], [L_y1], [L_z1]
    e1x, e1y, e1z = [0], [0], [0]
    e2x, e2y, e2z = [0], [0], [0]
    e3x, e3y, e3z = [0], [0], [0]
    for _ in range(steps):
        c1 = Lx_dot(L_y[-1], L_z[-1])
        d1 = Ly_dot(L_x[-1], L_z[-1])
        e1 = Lz_dot(L_x[-1], L_y[-1])
        f1 = e1x_dot()
        g1 = e1y_dot(L_z[-1])
        h1 = e1z_dot(L_y[-1])
        i1 = e2x_dot(L_z[-1])
        j1 = e2y_dot()
        k1 = e2z_dot(L_x[-1])
        l1 = e3x_dot(L_y[-1])
        m1 = e3y_dot(L_x[-1])
        n1 = e3z_dot()

        c2 = Lx_dot(L_y[-1] + (h * d1 / 2), L_z[-1] + (h * e1 / 2))
        d2 = Ly_dot(L_x[-1] + (h * c1 / 2), L_z[-1] + (h * e1 / 2))
        e2 = Lz_dot(L_x[-1] + (h * c1 / 2), L_y[-1] + (h * d1 / 2))
        f2 = e1x_dot()
        g2 = e1y_dot(L_z[-1] + (h * e1 / 2))
        h2 = e1z_dot(L_y[-1] + (h * d1 / 2))
        i2 = e2x_dot(L_z[-1] + (h * e1 / 2))
        j2 = e2y_dot()
        k2 = e2z_dot(L_x[-1] + (h * c1 / 2))
        l2 = e3x_dot(L_y[-1] + (h * d1 / 2))
        m2 = e3y_dot(L_x[-1] + (h * c1 / 2))
        n2 = e3z_dot()

        c3 = Lx_dot(L_y[-1] + (h * d2 / 2), L_z[-1] + (h * e2 / 2))
        d3 = Ly_dot(L_x[-1] + (h * c2 / 2), L_z[-1] + (h * e2 / 2))
        e3 = Lz_dot(L_x[-1] + (h * c2 / 2), L_y[-1] + (h * d2 / 2))
        f3 = e1x_dot()
        g3 = e1y_dot(L_z[-1] + (h * e2 / 2))
        h3 = e1z_dot(L_y[-1] + (h * d2 / 2))
        i3 = e2x_dot(L_z[-1] + (h * e2 / 2))
        j3 = e2y_dot()
        k3 = e2z_dot(L_x[-1] + (h * c2 / 2))
        l3 = e3x_dot(L_y[-1] + (h * d2 / 2))
        m3 = e3y_dot(L_x[-1] + (h * c2 / 2))
        n3 = e3z_dot()

        c4 = Lx_dot(L_y[-1] + (h * d3), L_z[-1] + (h * e3))
        d4 = Ly_dot(L_x[-1] + (h * c3), L_z[-1] + (h * e3))
        e4 = Lz_dot(L_x[-1] + (h * c3), L_y[-1] + (h * d3))
        f4 = e1x_dot()
        g4 = e1y_dot(L_z[-1] + (h * e3))
        h4 = e1z_dot(L_y[-1] + (h * d3))
        i4 = e2x_dot(L_z[-1] + (h * e3))
        j4 = e2y_dot()
        k4 = e2z_dot(L_x[-1] + (h * c3))
        l4 = e3x_dot(L_y[-1] + (h * d3))
        m4 = e3y_dot(L_x[-1] + (h * c3))
        n4 = e3z_dot()

        L_x.append(L_x[-1] + ((h / 6) * (c1 + (2 * c2) + (2 * c3) + c4)))
        L_y.append(L_y[-1] + ((h / 6) * (d1 + (2 * d2) + (2 * d3) + d4)))
        L_z.append(L_z[-1] + ((h / 6) * (e1 + (2 * e2) + (2 * e3) + e4)))
        e1x.append(e1x[-1] + ((h / 6) * (f1 + (2 * f2) + (2 * f3) + f4)))
        e1y.append(e1y[-1] + ((h / 6) * (g1 + (2 * g2) + (2 * g3) + g4)))
        e1z.append(e1z[-1] + ((h / 6) * (h1 + (2 * h2) + (2 * h3) + h4)))
        e2x.append(e2x[-1] + ((h / 6) * (i1 + (2 * i2) + (2 * i3) + i4)))
        e2y.append(e2y[-1] + ((h / 6) * (j1 + (2 * j2) + (2 * j3) + j4)))
        e2z.append(e2z[-1] + ((h / 6) * (k1 + (2 * k2) + (2 * k3) + k4)))
        e3x.append(e3x[-1] + ((h / 6) * (l1 + (2 * l2) + (2 * l3) + l4)))
        e3y.append(e3y[-1] + ((h / 6) * (m1 + (2 * m2) + (2 * m3) + m4)))
        e3z.append(e3z[-1] + ((h / 6) * (n1 + (2 * n2) + (2 * n3) + n4)))

    return np.array(L_x), np.array(L_y), np.array(L_z)

A = 1
B = 2
C = 3

def main():
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})

    L_x1_list = np.arange(-1, 1.5, 0.5)
    L_y1_list = np.arange(-2, 2, 1)
    L_z1_list = np.arange(-3, 3, 1.5)

    # Plot the sphere for visual reference
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    sphere_x = 0.98 * np.outer(np.cos(u), np.sin(v))
    sphere_y = 0.98 * np.outer(np.sin(u), np.sin(v))
    sphere_z = 0.98 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(
        sphere_x, sphere_y, sphere_z, color='r', alpha=0.2, edgecolor='none'
    )

    for L_x1 in L_x1_list:
        for L_y1 in L_y1_list:
            for L_z1 in L_z1_list:
                L_x, L_y, L_z = runge_kutta(L_x1, L_y1, L_z1)
                L = np.sqrt((L_x**2) + (L_y**2) + (L_z**2))
                L_x = L_x / L
                L_y = L_y / L
                L_z = L_z / L
                ax.plot(L_x, L_y, L_z, color="black")


    ax.set_xlabel("Lx")
    ax.set_ylabel("Ly")
    ax.set_zlabel("Lz")
    ax.set_title("Path of the L-point")
    ax.legend()

    plt.show()

if __name__ == "__main__":
    main()
