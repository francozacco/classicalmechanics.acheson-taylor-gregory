import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np


def main():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    X = np.arange(0, 6, 0.25)
    Y = np.arange(0, 6, 0.25)
    X, Y = np.meshgrid(X, Y)
    Z = (Y * np.e**(-Y)) * (X * np.e**(-X))

    # Plot the surface.
    surf = ax.plot_surface(
        X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False
    )

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


if __name__ == "__main__":
    main()
