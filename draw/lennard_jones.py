import numpy as np
import matplotlib.pyplot as plt
import os

path = os.path.dirname(__file__)

L = 10
r_cut = L/2

lj = lambda r : 4 * (r**-12 - r**-6)
lj_cut = lambda r, r_cut: lj(r) * np.array(r < r_cut, dtype=float)

fig, ax = plt.subplots()

rlin = np.linspace(1e-1, L, 1000)
rlin_ins = np.linspace(4.9, 5.1, 1000)
vert = np.linspace(-1, 1, 1000)

ax.set_xlim(left = np.min(rlin), right = np.max(rlin))
ax.set_ylim(bottom = np.min(vert), top = np.max(vert))

ax.plot(rlin, lj_cut(rlin, r_cut), ls = "solid", color = "tab:blue", label = r"$\phi_c^*$")
ax.plot(rlin, lj(rlin), ls = "dashed", color = "tab:orange", label = r"$\phi^*$")

ax.plot(r_cut * np.ones_like(vert), vert, ls = "dashed", color = "tab:gray", zorder=-1)

ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
ax.xaxis.set_minor_locator(plt.MultipleLocator(.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(.5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(.25))

ax.set_xlabel(r"$r_{ij}^*$")
ax.set_ylabel(r"$\phi_{(c)}^*$")

ax_ins = ax.inset_axes(
    [0.65, 0.1, 0.3, 0.3],
    xlim = (4.9, 5.1), ylim = (-5e-4, 0)
)

ax_ins.plot(rlin_ins, lj_cut(rlin_ins, r_cut), ls = "solid", color = "tab:blue", label = r"$\phi_c^*$")
ax_ins.plot(rlin_ins, lj(rlin_ins), ls = "dashed", color = "tab:orange", label = r"$\phi^*$")

ax_ins.plot(r_cut * np.ones_like(vert), vert, ls = "dashed", color = "tab:gray", zorder=-1)

ax.indicate_inset_zoom(ax_ins, edgecolor = "black")

ax.legend(loc = "best")

fig.savefig(path + "/lennard_jones.pdf", dpi = 300, bbox_inches = "tight")
