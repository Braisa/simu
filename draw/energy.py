import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

E = -680
logging_interval = 100

name = input("Name? ")

energy = pd.read_csv(path + name + ".txt", sep=" ")["E"]

fig, ax = plt.subplots()

steps = logging_interval * len(energy)
steps_lin = logging_interval * (1 + np.arange(len(energy)))

ax.plot(steps_lin, energy, ls = "solid", color = "tab:orange")
#ax.plot(steps_lin, np.mean(energy) * np.ones_like(steps_lin), ls = "dashed", color = "tab:orange", label="All steps mean")
#ax.plot(steps_lin, np.mean(energy[20:]) * np.ones_like(steps_lin), ls = "dashed", color = "tab:purple", label="Last steps mean")

ax.set_xlim(left=np.min(steps_lin), right=np.max(steps_lin))

ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))
ax.yaxis.set_major_locator(plt.MultipleLocator(.02))
ax.yaxis.set_minor_locator(plt.MultipleLocator(.01))

ax.set_xlabel("Steps")
ax.set_ylabel("E")

ax.legend(loc="best")

fig.savefig(path + name + ".pdf", dpi=300, bbox_inches="tight")
