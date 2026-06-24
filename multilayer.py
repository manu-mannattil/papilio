"""Multilayer interference calculations.

- Variable conventions and equations correspond to those in
  Peatross and Ware (2015) [1].
- Requires the Python package colour-science for CIE color matching
  functions.

TODO
----

1. There's a potential branch cut issue with arcsin().

[1] https://optics.byu.edu/docs/OpticsBook.pdf
"""

import numpy as np
import colour

nm = 1e-9
um = 1e-6

def reflectance(lam, n_list, d_list, theta_0=0, n_0=1, n_sub=1):
    """Compute the multilayer reflectance for unpolarized light.

    Parameters
    ----------
    lam : float
        Wavelength (in meters).
    n_list : array-like
        Refractive indices of layers.
    d_list : array-like
        Thickness of layers (in meters)
    theta_0 : float
        Initial incidence angle, default = 0.
    n_0 : float
        Refractive index of initial medium (i.e., outside the
        multilayer), default = 1.
    n_sub : float
        Refractive index of substrate, default = 1.

    Returns
    -------
    r : float
        Reflectance
    """
    r_p = reflectance_p(lam, n_list, d_list, theta_0, n_0, n_sub)
    r_s = reflectance_s(lam, n_list, d_list, theta_0, n_0, n_sub)
    return 0.5 * (r_p + r_s)

def reflectance_p(lam, n_list, d_list, theta_0=0, n_0=1, n_sub=1):
    """Compute the multilayer reflectance for p-polarized light."""
    # Total number of layers.
    N = len(n_list)

    # Compute the incidence angles in each layer using Snell's law.
    theta_list = np.empty(N, dtype=complex)
    theta_list[0] = np.arcsin(n_0 * np.sin(theta_0) / n_list[0])
    for i in range(1, N):
        theta_list[i] = np.arcsin(n_list[i - 1] * np.sin(theta_list[i - 1]) / n_list[i])

    # Exit angle (into the substrate).
    theta_sub = np.arcsin(n_list[N - 1] * np.sin(theta_list[N - 1]) / n_sub)

    # Phase shifts in each layer (Eq. 4.49).
    beta_list = np.empty(N, dtype=complex)
    for i in range(N):
        beta_list[i] = (2 * np.pi * n_list[i] / lam) * d_list[i] * np.cos(theta_list[i])

    # Equation 4.55.
    M = np.eye(2, dtype=complex)
    for i in range(N):
        theta, beta, n = theta_list[i], beta_list[i], n_list[i]
        M @= np.array([[np.cos(beta), -1j * np.sin(beta) * np.cos(theta) / n],
                       [-1j * n * np.sin(beta) / np.cos(theta), np.cos(beta)]])

    # Equation 4.58.
    A = np.eye(2, dtype=complex)
    A /= 2 * n_0 * np.cos(theta_0)
    A @= np.array([[n_0, np.cos(theta_0)], [n_0, -np.cos(theta_0)]])
    A @= M
    A @= np.array([[np.cos(theta_sub), 0], [n_sub, 0]])

    # Reflectance (Eq. 4.60).
    r = A[1][0] / A[0][0]
    return np.abs(r)**2

def reflectance_s(lam, n_list, d_list, theta_0=0, n_0=1, n_sub=1):
    """Compute the multilayer reflectance for s-polarized light."""
    # Total number of layers.
    N = len(n_list)

    # Compute the incidence angles in each layer using Snell's law.
    theta_list = np.empty(N, dtype=complex)
    theta_list[0] = np.arcsin(n_0 * np.sin(theta_0) / n_list[0])
    for i in range(1, N):
        theta_list[i] = np.arcsin(n_list[i - 1] * np.sin(theta_list[i - 1]) / n_list[i])

    # Exit angle (into the substrate).
    theta_sub = np.arcsin(n_list[N - 1] * np.sin(theta_list[N - 1]) / n_sub)

    # Phase shifts in each layer (Eq. 4.49).
    beta_list = np.empty(N, dtype=complex)
    for i in range(N):
        beta_list[i] = (2 * np.pi * n_list[i] / lam) * d_list[i] * np.cos(theta_list[i])

    # Equation 4.63.
    M = np.eye(2, dtype=complex)
    for i in range(N):
        theta, beta, n = theta_list[i], beta_list[i], n_list[i]
        M @= np.array([[np.cos(beta), -1j * np.sin(beta) / (n * np.cos(theta))],
                       [-1j * n * np.cos(theta) * np.sin(beta), np.cos(beta)]])

    # Equation 4.62.
    A = np.eye(2, dtype=complex)
    A /= 2 * n_0 * np.cos(theta_0)
    A @= np.array([[n_0 * np.cos(theta_0), 1], [n_0 * np.cos(theta_0), -1]])
    A @= M
    A @= np.array([[1, 0], [n_sub * np.cos(theta_sub), 0]])

    # Reflectance (Eq. 4.65).
    r = A[1][0] / A[0][0]
    return np.abs(r)**2

def reflectance_to_hex(lam_list,
                       r_list,
                       illuminant="D65",
                       observer="CIE 1931 2 Degree Standard Observer",
                       clip=True,
                       ):
    """
    Convert a reflectance-vs-wavelength curve to a displayable sRGB hex color
    using CIE 1931 color matching functions.

    Written by GPT 5.5.

    Parameters
    ----------
    lam_list : array-like
        Wavelength samples in meters.
    r_list : array-like
        Reflectance values at those wavelengths. Usually in [0, 1].
    illuminant : str
        Standard illuminant name, default "D65".
    observer : str
        CIE observer/color matching function set.
    clip : bool
        If True, clips out-of-gamut sRGB values to [0, 1].

    Returns
    -------
    color : str
        Hex color, e.g. "#8FA34C".
    """
    # Convert wavelengths to nm.
    lam_list = np.asarray(lam_list, dtype=float) / nm
    r_list = np.asarray(r_list, dtype=float)

    if np.any(np.diff(lam_list) <= 0):
        raise ValueError("lam_list must be strictly increasing.")

    # Build spectral reflectance distribution.
    sd = colour.SpectralDistribution(dict(zip(lam_list, r_list)), name="Sample reflectance", )

    # Use a 1 nm integration grid over the supplied wavelength range.
    shape = colour.SpectralShape(int(np.ceil(lam_list.min())), int(np.floor(lam_list.max())), 1, )

    sd = sd.copy().align(shape)

    cmfs = colour.MSDS_CMFS[observer].copy().align(shape)
    illuminant_sd = colour.SDS_ILLUMINANTS[illuminant].copy().align(shape)

    # Integrate reflectance * illuminant * CIE 1931 CMFs.
    # colour.sd_to_XYZ returns XYZ scaled so Y=100 for a perfect diffuser.
    XYZ = colour.sd_to_XYZ(sd, cmfs=cmfs, illuminant=illuminant_sd)

    # Convert XYZ [0, 100] to XYZ [0, 1].
    XYZ = XYZ / 100.0

    # Convert to display sRGB.
    rgb = colour.XYZ_to_sRGB(XYZ)

    if clip:
        rgb = np.clip(rgb, 0.0, 1.0)

    rgb8 = np.round(rgb * 255).astype(np.uint8)

    return "#{:02X}{:02X}{:02X}".format(*rgb8)
