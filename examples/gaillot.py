"""P. blumei reflectance spectra from Gaillot et al., Phys. Rev. E, 2008."""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import *

n_air = 1
n_sub = 1.52 # substrate refractive index (glass)
F = 70/170 # fill factor (Gaillot's estimate)

# Thickness list.
d_list = np.concatenate(([95*nm, 85*nm]*8, [95*nm]))

# Define wavelength region of interest.
M = 500
lam_list = np.linspace(400, 1200, M) * nm

# Incidence angles.
theta_0_list = np.linspace(0, 60, 7)

def n_chitin(lam):
    # Complex refractive index of chitin.
    # These values are a bit different from Gaillot et al.
    n = 1.7 + (1.5 - 1.7)/(800 - 400)*(lam - 400)
    k = 0.06 + (0.03 - 0.06)/(800 - 400)*(lam - 400)

    return n + 1j*k

def n_eff(lam):
    # Effective index of intermediate air layer.
    return (1 - F)*n_air + F*n_chitin(lam)

def _reflectance(theta_0):
    """Find the reflectance for a given incident angle."""
    r_list = np.empty(M)
    for i, lam in enumerate(lam_list):
        n1, n2 = n_chitin(lam), n_eff(lam)
        n_list = np.concatenate(([n1, n2]*8, [n1]))
        r_list[i] = reflectance(lam, n_list, d_list, theta_0, n_air, n_sub)

    return r_list

r_list_avg = np.zeros(M)
for theta_0 in theta_0_list:
    r_list = _reflectance(theta_0*np.pi/180)
    r_list_avg += r_list
    c = reflectance_to_hex(lam_list, r_list)
    plt.plot(lam_list / nm, r_list*100, color=c, label=f"{theta_0} deg")

r_list_avg /= len(theta_0_list)
c_avg = reflectance_to_hex(lam_list, r_list_avg)
plt.plot([], [], color=c_avg, label="average color", linewidth=10)

plt.xlim(400, 1200)
plt.title("p. blumei reflectance")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.legend()
plt.show()
