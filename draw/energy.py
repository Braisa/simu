import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

E = -680

energy = pd.read_csv(path + "energy.txt", sep=" ")["E"]

fig, ax = plt.subplots()

steps = 1 + np.arange(len(energy))

ax.plot(steps, energy, ls = "solid", color = "tab:blue")
ax.plot(steps, np.mean(energy) * np.ones_like(steps), ls = "dashed", color = "tab:orange", label="All steps mean")
ax.plot(steps, np.mean(energy[20:]) * np.ones_like(steps), ls = "dashed", color = "tab:purple", label="Last steps mean")

ax.set_xlim(left=np.min(steps), right=np.max(steps))

ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2))
ax.yaxis.set_major_locator(plt.MultipleLocator(.02))
ax.yaxis.set_minor_locator(plt.MultipleLocator(.01))

ax.set_xlabel("Steps")
ax.set_ylabel("E")

ax.legend(loc="best")

fig.savefig(path + "energy.pdf", dpi=300, bbox_inches="tight")
