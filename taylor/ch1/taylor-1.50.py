import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math


def main():
    plt.figure()
    ax = plt.axes()
    R = 5
    phi0 = 20
    omega_0 = math.sqrt(9.81 / R)
    phi, omega, t = improved_euler(phi0=phi0, omega0=0, R=R)
    plot(ax, t, phi, "numerical")
    plot(ax, t, [phi0 * math.cos(omega_0 * i) for i in t], "approx")
    plt.axvline(x=4.5, linewidth=1, color="grey", linestyle="--")
    plt.axvline(x=9, linewidth=1, color="grey", linestyle="--")
    plt.show()
    # plt.savefig("taylor-1.51.png")


def fphi(phi, omega):
    return omega * 180 / math.pi


def fomega(phi, omega, R):
    return - (9.81 / R) * math.sin(math.radians(phi))


def improved_euler(phi0, omega0, R, h=0.01, it=1000):
    phi, omega = [phi0], [omega0]
    t = [0]

    for _ in range(it):
        c1 = h * fphi(phi[-1], omega[-1])
        d1 = h * fomega(phi[-1], omega[-1], R)

        c2 = h * fphi(phi[-1] + c1, omega[-1] + d1)
        d2 = h * fomega(phi[-1] + c1, omega[-1] + d1, R)

        phi.append(phi[-1] + (1/2 * (c1 + c2)))
        omega.append(omega[-1] + (1/2 * (d1 + d2)))
        t.append(t[-1] + h)

    return phi, omega, t


def plot(ax, x, y, label):
    ax.plot(x, y, label=label)
    ax.legend()
    ax.set_xlabel("t")
    ax.set_ylabel("phi")


if __name__ == "__main__":
    main()
