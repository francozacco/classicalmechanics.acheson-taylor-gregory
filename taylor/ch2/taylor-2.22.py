import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import numpy as np


def main():
    plt.figure()
    ax = plt.axes()

    R_list = []
    theta_list = np.arange(0.1, 0.8, 0.02)
    for theta in theta_list:
        f_mid, R_mid = binary_search(theta=theta)
        print(f"theta: {theta} f_mid: {f_mid} R_mid: {R_mid}")
        R_list.append(R_mid)
    plot(ax, theta_list, R_list, "R")
    plt.axvline(
        x=theta_list[np.argmax(R_list)],
        linewidth=1,
        color="grey",
        linestyle="--",
        label=f'theta: {round(theta_list[np.argmax(R_list)], 2)}\n R: {round(max(R_list), 2)}'
    )
    ax.legend()

    plt.show()
    # plt.savefig("taylor-2.22.png")


def func(R, theta):
    try:
        return (((math.sin(theta) + 1) / math.cos(theta)) * R) + \
            math.log(1 - (R / math.cos(theta)))
    except:
        import ipdb; ipdb.set_trace()

def binary_search(theta):

    R_max = 0.01
    R_min = 0.70
    epsilon = 0.000001
    f_mid = 1000

    while abs(f_mid) > epsilon:
        R_mid = R_min + ((R_max - R_min)/2)
        f_mid = func(R_mid, theta)
        if f_mid > 0:
            R_max = R_mid
        elif f_mid < 0:
            R_min = R_mid

    return f_mid, R_mid


def plot(ax, x, y, label):
    ax.plot(x, y, label=label)
    ax.set_xlabel("theta")
    ax.set_ylabel("R")


if __name__ == "__main__":
    main()
