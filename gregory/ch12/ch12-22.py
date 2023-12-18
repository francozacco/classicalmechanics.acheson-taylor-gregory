import matplotlib
import matplotlib.pyplot as plt
import math
matplotlib.use('TkAgg')


def dtheta_dt_fn(v):
    return v


def dv_dt_fn(epsilon, p_Omega, t, theta):
    return (((1 / p_Omega) ** 2) + (epsilon * math.cos(t))) * math.sin(theta)


def runge_kutta(ax1, ax2, epsilon, p_Omega, theta_init, v_init):
    h = 0.001
    steps = 200000

    theta = [theta_init]
    v = [v_init]
    t = [0]
    for _ in range(steps):
        c1 = dtheta_dt_fn(v[-1])
        d1 = dv_dt_fn(epsilon, p_Omega, t[-1], theta[-1])

        c2 = dtheta_dt_fn(v[-1] + (h * d1 / 2))
        d2 = dv_dt_fn(epsilon, p_Omega, t[-1] + (h / 2), theta[-1] + (h * c1 / 2))

        c3 = dtheta_dt_fn(v[-1] + (h * d2 / 2))
        d3 = dv_dt_fn(epsilon, p_Omega, t[-1] + (h / 2), theta[-1] + (h * c2 / 2))

        c4 = dtheta_dt_fn(v[-1] + (h * d3))
        d4 = dv_dt_fn(epsilon, p_Omega, t[-1] + (h / 2), theta[-1] + (h * c3))

        theta.append(theta[-1] + ((h / 6) * (c1 + (2 * c2) + (2 * c3) + c4)))
        v.append(v[-1] + ((h / 6) * (d1 + (2 * d2) + (2 * d3) + d4)))
        t.append(t[-1] + h)
    ax1.plot(t, v)
    ax2.plot(t, theta)


def main():
    _, ax1 = plt.subplots(figsize=(10, 10))
    _, ax2 = plt.subplots(figsize=(12, 7))
    legends = []
    epsilon = 0.3
    p_Omega = [1, 2, 5, 7, 10]
    theta_init_list = [0.1, 0.1, 0.1, 0.1, 0.1]
    v_init_list = [0, 0, 0, 0, 0]

    for p_O, theta_init, v_init in zip(p_Omega, theta_init_list, v_init_list):
        runge_kutta(ax1, ax2, epsilon, p_O, theta_init, v_init)
        legends.append(f"e:{epsilon} p/Omega:{p_O} theta_init:{theta_init}")

    ax1.legend(labels=legends)
    ax1.set_xlabel("t")
    ax1.set_ylabel("v")
    ax1.grid(True)

    ax2.legend(labels=legends)
    ax2.set_xlabel("t")
    ax2.set_ylabel("theta")
    ax2.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
