import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    plt.figure()
    ax = plt.axes()

    omega = 30
    epsilon = 0.01
    F0 = 300
    t_range = np.arange(0, 40, 0.1)

    x = []
    for t in t_range:
        c = F0/((omega**2)*epsilon*(1+(epsilon/2)))
        x.append(c*np.sin(omega*t*(1+(epsilon/2)))*np.sin(omega*t*epsilon/2))
    textstr = (
        r'$\Omega = $' + f"{omega}\n" + 
        r'$\epsilon = $' + f"{epsilon}\n" +
        r'$F_0 = $' + f"{F0}\n"
    )
    plot(ax, t_range, x, textstr)
    plt.show()
    # plt.savefig("gregory-5.10.png")

def plot(ax, x, y, textstr):
    ax.plot(x, y)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.text(
        0.05,
        0.95,
        textstr,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment='top'
    )

if __name__ == "__main__":
    main()