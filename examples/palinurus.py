"""P. palinurus reflectance curves.

The idea is that the scales are so curved that rays undergo a double
reflection.

Ref.: Vukusic et al., Appl. Opt., 2001.
"""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import *

fig, ax = plt.subplots()

lam_obs, r_obs = np.loadtxt("data/palinurus.dat", unpack=True)
ax.plot(lam_obs, r_obs, "o", zorder=100)

N_layers = 10 # number of layers
n_air = 1.0 # index of air
n_sub = 1.52 # substrate refractive index (glass)
n_chitin = 1.58 # guess
F = 0.5 # fill factor (guess)
k_min, k_max = 0.06, 0.25 # extinction coefficient (Vukusic's estimate)

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

# Double reflection.
# The intensity has to be squared.
r_45 = np.array([reflectance(lam, n_list, d_list, np.pi/4, n_air, n_sub) for lam in lam_list])
r_45 = r_45**2
c_45 = reflectance_to_hex(lam_list, r_45)
ax.plot(lam_list, r_45*100, color=c_45, label="45 deg")

# Normal incidence.
r_90 = np.array([reflectance(lam, n_list, d_list, 0, n_air, n_sub) for lam in lam_list])
c_90 = reflectance_to_hex(lam_list, r_90)
ax.plot(lam_list, r_90*100, color=c_90, label="0 deg")

r_both = r_90 + r_45
c_both = reflectance_to_hex(lam_list, r_both)
ax.plot(lam_list, 100*r_both, color=c_both, label="0 + 45 deg")
ax.plot([], [], color=c_both, linewidth=10, label="perceived color")

# Fig. 8 of Vukusic et al.
ax.set_xlim(325, 725)
ax.set_title(r"P. palinurus reflectance")
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("reflectance (%)")
ax.legend()
plt.show()
