import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')


def main():
    eps = 0.01

    rho_init = 1
    sigma_init = 0
    theta_init = 0

    h = 0.0001
    steps = 500000

    rho = [rho_init]
    sigma = [sigma_init]
    theta = [theta_init]
    t = 0
    for _ in range(steps):
        c1 = rho_fn(sigma[-1])
        d1 = sigma_fn(eps, t, rho[-1], sigma[-1])
        f1 = theta_fn(eps, t, rho[-1])

        c2 = rho_fn(sigma[-1] + (h * c1 / 2))
        d2 = sigma_fn(eps, t + (h / 2), rho[-1] + (h * d1 / 2), sigma[-1] + (h * d1 / 2))
        f2 = theta_fn(eps, t + (h / 2), rho[-1] + (h * f1 / 2))

        c3 = rho_fn(sigma[-1] + (h * c2 / 2))
        d3 = sigma_fn(eps, t + (h / 2), rho[-1] + (h * d2 / 2), sigma[-1] + (h * d2 / 2))
        f3 = theta_fn(eps, t + (h / 2), rho[-1] + (h * f2 / 2))

        c4 = rho_fn(sigma[-1] + (c3 * h))
        d4 = sigma_fn(eps, t + h, rho[-1] + (h * d2), sigma[-1] + (h * d2))
        f4 = theta_fn(eps, t + h, rho[-1] + (h * f3))

        rho.append(rho[-1] + (h / 6 * (c1 + (2 * c2) + (2 * c3) + c4)))
        sigma.append(sigma[-1] + (h / 6 * (d1 + (2 * d2) + (2 * d3) + d4)))
        theta.append(theta[-1] + (h / 6 * (f1 + (2 * f2) + (2 * f3) + f4)))
        t += h

    plot(theta, rho)
    plt.show()
    # plt.savefig("gregory-7.28.png")


def rho_fn(sigma):
    return sigma


def sigma_fn(eps, t, rho, sigma):
    return (math.e**(-2 * eps * t) / rho**3) - (1 / rho**2) - (eps * sigma)


def theta_fn(eps, t, rho):
    return math.e**(-2 * eps * t) / rho**2


def plot(theta, r):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, r)
    ax.set_rmax(1)
    ax.grid(True)
    # ax.set_title("Mercury Perihelion Precession")


if __name__ == "__main__":
    main()
