"""P. blumei reflectance spectra from Tada et al., Appl. Opt., 1998."""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import *

n_air = 1
n_sub = 1.52 # substrate refractive index (glass)
n_chitin = 1.58 # see paragraph after Fig. 13.
F = 0.5 # fill factor (Gaillot's estimate)
n_eff = (1 - F)*n_air + F*n_chitin

# Thickness and index list.
d_list = np.array([95, 85]*10)
n_list = np.array([n_chitin, n_eff]*10)

# Define wavelength region of interest.
M = 500
lam_list = np.linspace(300, 900, M)

# Incidence angles.
theta_0_list = np.linspace(0, 60, 7)

def _reflectance(theta_0):
    # Find the reflectance for a given incidence angle
    r_list = np.empty(M)
    for i, lam in enumerate(lam_list):
        r_list[i] = reflectance(lam, n_list, d_list, theta_0, n_air, n_sub)

    return r_list

r_list_avg = np.zeros(M)
for theta_0 in theta_0_list:
    r_list = _reflectance(theta_0*np.pi/180)
    r_list_avg += r_list
    c = reflectance_to_hex(lam_list, r_list)
    plt.plot(lam_list, r_list*100, color=c, label=f"{theta_0} deg")

r_list_avg /= len(theta_0_list)
c_avg = reflectance_to_hex(lam_list, r_list_avg)
plt.plot([], [], color=c_avg, label="average color", linewidth=10)

# Fig. 13 of Tada et al.
plt.xlim(300, 900)
plt.title("p. blumei reflectance")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.legend()
plt.show()
