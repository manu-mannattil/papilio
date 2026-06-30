"""P. ulysses reflectance curves.

Ref.: Vukusic et al., Appl. Opt., 2001.
"""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import *
from scipy.optimize import least_squares, curve_fit

fig, ax = plt.subplots()

lam_obs, r_obs = np.loadtxt("data/ulysses.dat", unpack=True)
ax.plot(lam_obs, r_obs, "o", zorder=100)

N_layers = 10 # number of layers
n_air = 1.0 # index of air
n_sub = 1.52 # substrate refractive index (glass)
n_chitin = 1.58 # guess
F = 0.1 # fill factor (guess)
k_min, k_max = 0.02, 0.3 # extinction coefficient (guess)

# Thickness and index list.
d_list = np.array([100, 90]*N_layers)
k_list = np.linspace(k_min, k_max, N_layers)
n_chitin = np.array(n_chitin + 1j*k_list)
n_eff = (1 - F)*n_air + F*n_chitin
n_list = np.vstack([n_chitin, n_eff])
n_list = n_list.T.ravel()

# Define wavelength region of interest.
M = 100
lam_list = np.linspace(325, 725, M)

# Compute the average reflectance by averaging over a range of angles.
def _reflectance(lam_list, theta_max=60*np.pi/180, steps=7):
    r_list_avg = np.zeros(M)
    for theta_0 in np.linspace(0, theta_max, steps):
        r_list = np.array([reflectance(lam, n_list, d_list, theta_0, n_air, n_sub) for lam in lam_list])
        r_list_avg += r_list

    return r_list_avg/steps

r_avg = _reflectance(lam_list)
color = reflectance_to_hex(lam_list, r_avg)
ax.plot(lam_list, 100*r_avg, color=color)
ax.plot([], [], color=color, linewidth=10, label="perceived color")

# Fig. 13 of Tada et al.
ax.set_xlim(325, 725)
ax.set_title(r"P. ulysses reflectance")
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("reflectance (%)")
ax.legend()
plt.show()
