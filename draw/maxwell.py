import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

rva = pd.read_csv(path + "eq_readable.txt", sep=" ")

rva["v2"] = rva["vx"]**2 + rva["vy"]**2 + rva["vz"]**2

names = ("vx", "vy", "vz", "v2")
xlabels = (r"$v_x$", r"$v_y$", r"$v_z$", r"$v^2$")
ylabels = (r"$f(v_x)$", r"$f(v_y)$", r"$f(v_z)$", r"$f(v)$")
colors = ("red", "green", "blue", "orange")

for (name, xlabel, ylabel, color) in zip(names, xlabels, ylabels, colors):
    fig, ax = plt.subplots()

    if not name == "v2":
        ax.set_xlim(left=-4, right=+4)
        b = 800
    else:
        ax.set_xlim(left=0, right=25)
        b = 2500

    ax.hist(rva[name], bins=b, density=True, histtype="step", color=f"tab:{color}")
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(path + f"maxwell_{name}.pdf", dpi=300, bbox_inches="tight")
