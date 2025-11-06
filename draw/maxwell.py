import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

rva = pd.read_csv(path + "big_rva.txt", sep=" ")

rva["v2"] = rva["vx"]**2 + rva["vy"]**2 + rva["vz"]**2

names = ("vx", "vy", "vz", "v2")
xlabels = (r"$x$", r"$y$", r"$z$", r"$r^2$")
ylabels = (r"$v_x$", r"$v_y$", r"$v_z$", r"$v^2$")
colors = ("red", "green", "blue", "orange")

for (name, xlabel, ylabel, color) in zip(names, xlabels, ylabels, colors):
    fig, ax = plt.subplots()

    plt.hist(rva[name], bins=50, histtype="step", color=f"tab:{color}")
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if not name == "v2":
        ax.set_xlim(left=-4, right=+4)
    else:
        ax.set_xlim(left=0, right=25)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(path + f"maxwell_{name}.pdf", dpi=300, bbox_inches="tight")
