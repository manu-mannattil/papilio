"""Broadband reflector example from Jun Zhang's notes [1].

[1] https://wp.optics.arizona.edu/milster/wp-content/uploads/sites/48/2016/06/Thin-film-calculator-from-Dissertation_JunZhang_080110_optimized.pdf
"""

from multilayer import *
import numpy as np
import matplotlib.pyplot as plt

# Incidence angle.
theta_0 = 0

# Refractive indices.
n_air = 1.0
n_sub = 1.52 # glass
nH = 2.35
nL = 1.35
n_list = np.array([nH, nL, nH, nL, nH, nL, nH, nL, nH, nL, nH, nL,
                   nH, nL, nH, nL, nH, nL, nH, nL, nH, nL, nH], dtype=complex)

# Define reference wavelength and quarter-wave layer thicknesses.
lam_f = 480 * nm
dH = lam_f / (4 * nH)
dL = lam_f / (4 * nL)

# Layer thicknesses.
# See Table I of Zhang's notes.
d_list = np.array([dH, dL, dH, dL, dH, dL, dH, dL, dH, dL, dH,
                   1.2 * dL, 1.4 * dH, 1.4 * dL, 1.4 * dH, 1.4 * dL, 1.4 * dH,
                   1.4 * dL, 1.4 * dH, 1.4 * dL, 1.4 * dH, 1.4 * dL, 1.4 * dH], dtype=float)

# Define visible wavelength region.
samples = 500
lam_list = np.linspace(350, 850, samples) * nm

r_list = np.empty(samples)
for i, lam in enumerate(lam_list):
    r_list[i] = reflectance_s(lam, n_list, d_list, theta_0, n_air, n_sub)

plt.figure()
plt.plot(lam_list / nm, r_list*100)
plt.title("A broadband reflector for the visible region")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectance (%)")
plt.axis([350, 850, 0, 105])
plt.grid(True)
plt.show()
