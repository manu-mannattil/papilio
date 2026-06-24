"""Multilayer interference calculations.

All equations correspond to Peatross and Ware (2015).
https://optics.byu.edu/docs/OpticsBook.pdf
"""

import numpy as np

def reflectance_s(lam, theta_0, n_0, n_sub, n_list, d_list):
    # Total number of layers.
    N = len(n_list)

    # Compute the incident angles in each layer using Snell's law.
    theta_list = np.empty(N, dtype=complex)
    theta_list[0] = np.arcsin(n_0*np.sin(theta_0)/n_list[0])
    for i in range(1, N):
        theta_list[i] = np.arcsin(n_list[i - 1]*np.sin(theta_list[i - 1])/n_list[i])

    # Exit angle (into the substrate).
    theta_sub = np.arcsin(n_list[N - 1]*np.sin(theta_list[N - 1])/n_sub)

    # Phase shifts in each layer.
    # Equation 4.49 of Peatross & Ware (2015).
    beta_list = np.empty(N, dtype=complex)
    for i in range(N):
        beta_list[i] = (2*np.pi*n_list[i]/lam) * d_list[i] * np.cos(theta_list[i])

    # Equation 4.63.
    M = np.eye(2, dtype=complex)
    for i in range(N):
        theta, beta, n = theta_list[i], beta_list[i], n_list[i]
        M @= np.array([[np.cos(beta), -1j*np.sin(beta)/(n*np.cos(theta))],
                      [-1j*n*np.cos(theta)*np.sin(beta), np.cos(beta)]])

    # Equation 4.62.
    A = np.eye(2, dtype=complex)
    A /= 2*n_0*np.cos(theta_0)
    A @= np.array([[n_0*np.cos(theta_0), 1],
                   [n_0*np.cos(theta_0), -1]])
    A @= M
    A @= np.array([[1, 0],
                   [n_sub*np.cos(theta_sub), 0]])

    # Reflectance (Eq. 4.65).
    r = A[1][0]/A[0][0]
    return np.abs(r)**2
