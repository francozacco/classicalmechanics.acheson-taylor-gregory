import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')


def x_fn(v):
    return v


def v_fn(eps, x, v):
    return -eps * ((x**2) - 1) * v - x


def runge_kutta(ax1, ax2, eps=2, x_init=1, v_init=1):
    h = 0.001
    steps = 500000

    x = [x_init]
    v = [v_init]
    t = [0]
    for _ in range(steps):
        c1 = x_fn(v[-1])
        d1 = v_fn(eps, x[-1], v[-1])

        c2 = x_fn(v[-1] + (h * d1 / 2))
        d2 = v_fn(eps, x[-1] + (h * c1 / 2), v[-1] + (h * d1 / 2))

        c3 = x_fn(v[-1] + (h * d2 / 2))
        d3 = v_fn(eps, x[-1] + (h * c2 / 2), v[-1] + (h * d2 / 2))

        c4 = x_fn(v[-1] + (h * d3))
        d4 = v_fn(eps, x[-1] + (h * c3), v[-1] + (h * d3))

        x.append(x[-1] + (h / 6 * (c1 + (2 * c2) + (2 * c3) + c4)))
        v.append(v[-1] + (h / 6 * (d1 + (2 * d2) + (2 * d3) + d4)))
        t.append(t[-1] + h)
    ax1.plot(x, v)
    ax2.plot(t, x)


def main():
    _, ax1 = plt.subplots(figsize=(10, 10))
    _, ax2 = plt.subplots(figsize=(12, 7))
    legends = []
    eps_list = [2, 5, 10, 6]
    x_init_list = [1, 4, 1, 1]
    v_init_list = [1, 1, 10, 0]
    for eps, x_init, v_init in zip(eps_list, x_init_list, v_init_list):
        runge_kutta(ax1, ax2, eps, x_init, v_init)
        legends.append(f"eps:{eps} x_init:{x_init} v_init:{v_init}")
    ax1.legend(labels=legends)
    ax1.set_xlabel("x")
    ax1.set_ylabel("v")
    ax1.grid(True)
    ax1.set_title("Van der Pol's phase paths")

    ax2.legend(labels=legends)
    ax2.set_xlim(xmin=0, xmax=45)
    ax2.set_xlabel("t")
    ax2.set_ylabel("x")
    ax2.grid(True)
    ax2.set_title("Van der Pol's solution")

    plt.show()


if __name__ == "__main__":
    main()
