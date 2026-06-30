"""Layered material from Vukusic & Stavenega, J. R. Soc. Interface (2009)."""

import numpy as np
import matplotlib.pyplot as plt
from papilio import *

n_air = 1
n_sub = 1
n1 = 1.56 + 0.06j
n2 = 1

# Thickness and index list.
d_list = np.array([80, 100]*5)
n_list = np.array([n1, n2]*5, dtype=complex)

# Define wavelength region of interest.
M = 500
lam_list = np.linspace(350, 850, M)

# Incidence angle.
theta_0 = 0

r_list = np.empty(M)
for i, lam in enumerate(lam_list):
    r_list[i] = reflectance(lam, n_list, d_list, theta_0, n_air, n_sub)

# Fig. 11 of Vukusic & Stavenega.
fig, ax = plt.subplots()
ax.plot(lam_list, r_list*100, "C3")
ax.set_xlim(350, 850)
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("reflectance (%)")
plt.show()
