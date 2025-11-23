import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

rva = pd.read_csv(path + "eq_readable.txt", sep=" ")

steps = 500_000
logs = 5_000
spl = 500 # steps per log
N = 500
a = 2

names = ("rx", "ry", "rz")
ylabels = (r"$\lambda_x$", r"$\lambda_y$", r"$\lambda_z$")
colors = ("red", "green", "blue")

for (name, ylabel, color) in zip(names, ylabels, colors):

    order = np.zeros(logs)
    for l in range(logs):
        order[l] = 1/N * np.sum(np.cos(4*np.pi/a * rva[name][l*spl:(l+1)*spl]))
    
    fig, ax = plt.subplots()

    logs_lin = 100 * (1 + np.arange(logs))

    ax.plot(logs_lin, order, ls="solid", color=f"tab:{color}")
    ax.plot(logs_lin, 1/np.sqrt(N)*np.ones_like(logs_lin), ls="dashed", color="tab:gray")
    ax.plot(logs_lin, -1/np.sqrt(N)*np.ones_like(logs_lin), ls="dashed", color="tab:gray")

    ax.set_xlim(left=np.min(logs_lin), right=np.max(logs_lin))

    ax.set_xlabel("Steps")
    ax.set_ylabel(ylabel)

    ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))

    fig.tight_layout()
    fig.savefig(path + f"order_parameter_{name}.pdf", dpi=300, bbox_inches="tight")

    ax.plot()

