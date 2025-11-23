import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os
from scipy.optimize import curve_fit

path = os.path.dirname(__file__) + "/"

rva = pd.read_csv(path + "eq_readable.txt", sep=" ")

rva["v"] = np.sqrt(rva["vx"]**2 + rva["vy"]**2 + rva["vz"]**2)

steps = 500_000
logs = 5_000
spl = 500 # steps per log
width = 0.01

# For each log, we fit the f(v) distribution to the expected Maxwell distribution in order
# to find its corresponding temperature, from which we will be able to calculate the expected
# Boltzmann H-function value.

low, high = 0, 5
b = int((high-low)/width)
vlin = low + width*np.arange(b)

maxwell_v = lambda v, T : 4*np.pi*v**2 / np.sqrt(2*np.pi*T)**3 * np.exp(-(v**2)/(2*T))

H = np.zeros(logs)
H_maxwell = np.zeros(logs)

for l in range(logs):

    hist, _ = np.histogram(rva["v"][l*spl:(l+1)*spl], bins=b, range=(low,high), density=True)

    temp = curve_fit(maxwell_v, vlin, hist, p0=(1.))[0][0]

    hist_maxwell, _ = np.histogram(maxwell_v(vlin, temp), bins=b, range=(low,high), density=True)

    for (h, h_m) in zip(hist, hist_maxwell):
        if h != 0:
            H[l] += width * h*np.log(h)
        if h_m != 0:
            H_maxwell[l] += width * h_m*np.log(h_m)

fig, ax = plt.subplots()

logs_lin = 100 * (1 + np.arange(logs))

ax.plot(logs_lin, H, ls="solid", color="tab:orange")
ax.plot(logs_lin, H_maxwell, ls="solid", color="tab:gray")

ax.set_xlim(left=np.min(logs_lin), right=np.max(logs_lin))

ax.set_xlabel("Steps")
ax.set_ylabel(r"$H$")

ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))

fig.tight_layout()
fig.savefig(path + "boltzmann_H.pdf", dpi=300, bbox_inches="tight")

# Instead of using the v module distribution, we use each of the components of velocity

low, high = -4, 4
b = int((high-low)/width)
vilin = low + width*np.arange(b)

maxwell_v_i = lambda v_i, T : 1/np.sqrt(2*np.pi*T) * np.exp(-v_i**2/2/T)
comps = ("vx", "vy", "vz")

H = np.zeros(logs)
H_maxwell = np.zeros(logs)

for comp in comps:
    for l in range(logs):

        hist, _ = np.histogram(rva[comp][l*spl:(l+1)*spl], bins=b, range=(low,high), density=True)

        temp = curve_fit(maxwell_v_i, vilin, hist, p0=(1.))[0][0]

        hist_maxwell, _ = np.histogram(maxwell_v_i(vilin, temp), bins=b, range=(low,high), density=True)

        for (h, h_m) in zip(hist, hist_maxwell):
            if h != 0:
                H[l] += width/3 * h*np.log(h)
            if h_m != 0:
                H_maxwell[l] += width/3 * h_m*np.log(h_m)

fig, ax = plt.subplots()

logs_lin = 100 * (1 + np.arange(logs))

ax.plot(logs_lin, H, ls="solid", color="tab:orange")
ax.plot(logs_lin, H_maxwell, ls="solid", color="tab:gray")

ax.set_xlim(left=np.min(logs_lin), right=np.max(logs_lin))

ax.set_xlabel("Steps")
ax.set_ylabel(r"$H$")

ax.xaxis.set_major_locator(plt.MultipleLocator(steps//5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(steps//20))

fig.tight_layout()
fig.savefig(path + "boltzmann_H_comps.pdf", dpi=300, bbox_inches="tight")
