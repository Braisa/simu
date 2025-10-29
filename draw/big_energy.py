import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

E = -680
logging_interval = 100

name = input("Name? ")

big_energies = pd.read_csv(path + name + ".txt", sep=" ")

steps = logging_interval * len(big_energies["E"])
steps_lin = logging_interval * (1 + np.arange(len(big_energies["E"])))

fig, ax = plt.subplots()

ax.plot(steps_lin, big_energies["E"], ls="solid", color="tab:orange")
ax.plot(steps_lin, 1.001*E*np.ones_like(steps_lin), ls="dashed", color="tab:gray")
ax.plot(steps_lin, 0.999*E*np.ones_like(steps_lin), ls="dashed", color="tab:gray")

ax.set_xlim(left=np.min(steps_lin), right=np.max(steps_lin))

ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))
ax.yaxis.set_major_locator(plt.MultipleLocator(.0005*np.abs(E)))
ax.yaxis.set_minor_locator(plt.MultipleLocator(.0001*np.abs(E)))

ax.set_xlabel("Steps")
ax.set_ylabel("E")

fig.savefig(path + name + "_E.pdf", dpi=300, bbox_inches="tight")

fig, ax = plt.subplots()

ax.plot(steps_lin, big_energies["E"], ls="solid", color="tab:orange", label="E")
ax.plot(steps_lin, big_energies["V"], ls="solid", color="tab:purple", label="V")
ax.plot(steps_lin, big_energies["T"], ls="solid", color="tab:blue", label="T")

ax.set_xlim(left=np.min(steps_lin), right=np.max(steps_lin))

ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))

ax.set_xlabel("Steps")
ax.set_ylabel("Energies")

ax.legend(loc="best")

fig.savefig(path + name + ".pdf", dpi=300, bbox_inches="tight")
