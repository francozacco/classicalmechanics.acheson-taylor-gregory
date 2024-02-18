import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.use('TkAgg')


def plot_lambda_fn():
    _, ax1 = plt.subplots(figsize=(10, 10))

    _lambda = np.arange(-1, 5, 0.1)

    factor = 2

    fn_cosh = np.cosh(_lambda)
    fn_x = factor * _lambda

    ax1.plot(_lambda, fn_cosh)
    ax1.plot(_lambda, fn_x)

    ax1.legend(labels=[
        f"cosh(lambda) | b/a = {factor}",
        f"(b/a) * lambda | b/a = {factor}"
    ])
    ax1.set_xlabel("lambda")
    ax1.grid(True)

    plt.show()

def plot_extremals():
    _, ax1 = plt.subplots(figsize=(10, 10))

    x = np.arange(-3, 3, 0.1)

    C_1 = 1/0.5893
    C_2 = 1/2.1268

    y_1 = C_1 * np.cosh(x/C_1)
    y_2 = C_2 * np.cosh(x/C_2)

    ax1.plot(x, y_1)
    ax1.plot(x, y_2)

    ax1.legend(labels=[
        f"y1 | C_1 = 1/0.5893",
        f"y2 | C_2 = 1/2.1268"
    ])
    ax1.set_xlabel("x")
    ax1.set_ylim(0, 50)
    ax1.grid(True)

    plt.show()
def main():
    # plot_lambda_fn()
    plot_extremals()


if __name__ == "__main__":
    main()
