"""Reflections in a spherical cap."""

import numpy as np
import matplotlib.pyplot as plt

LOSS_RED = 0.75
LOSS_BLUE = 0.5

def z_and_R(theta):
    """Get complex representation of the initial ray and the rotation matrix."""
    z = -np.cos(np.pi/2 - theta) + 1j*np.sin(np.pi/2 - theta)
    R = np.cos(np.pi - 2*theta) + 1j*np.sin(np.pi - 2*theta)

    return z, R

def n_fold(n, theta, max_y=0.2):
    """Get piecewise rays for an n-fold reflection."""
    x, y = [], []

    z, R = z_and_R(theta)
    for i in range(n + 2):
        w = z*R**i
        x.append(w.real)
        y.append(w.imag)

    x, y = np.array(x), np.array(y)

    # Adjust end points so that y_max condition is satisfied.
    y[0] = max_y
    x[-1] = (max_y - y[-2])/(y[-1] - y[-2])*(x[-1] - x[-2]) + x[-2]
    y[-1] = max_y

    xx = np.array([x[:-1], x[1:]]).T
    yy = np.array([y[:-1], y[1:]]).T
    return xx, yy

# Plot the bowl.
fig, ax = plt.subplots()
theta = np.linspace(-np.pi, 0, 100)
ax.plot(np.cos(theta), np.sin(theta), "k", zorder=100)

# 45 < theta < 60.
x_list = np.linspace(np.sin(45*np.pi/180), np.sin(60*np.pi/180), 7)[1:-1]
dx = x_list[1] - x_list[0]
for theta in np.arcsin(x_list):
    if theta*180/np.pi > 54:
        n = 3
    else:
        n = 2
    xx, yy = n_fold(n, theta=theta)
    for n, (x, y) in enumerate(zip(xx, yy)):
        ax.plot(x, y, "C0", alpha=LOSS_BLUE**n)

# theta = 45, theta = 60.
xx, yy = n_fold(2, theta=np.pi/4)
for n, (x, y) in enumerate(zip(xx, yy)):
    ax.plot(x, y, "C3", alpha=LOSS_RED**n, zorder=50)
xx, yy = n_fold(3, theta=np.pi/3)
for n, (x, y) in enumerate(zip(xx, yy)):
    ax.plot(x, y, "C3", alpha=LOSS_RED**n, zorder=50)

# 30 <= theta < 45.
x_list = np.arange(np.sin(45*np.pi/180), np.sin(30*np.pi/180), -dx)[1:-1]
for theta in np.arcsin(x_list):
    xx, yy = n_fold(2, theta=theta)
    for n, (x, y) in enumerate(zip(xx, yy)):
        ax.plot(x, y, "C0", alpha=LOSS_BLUE**n)

# 60 < theta < 64 + 2/7
x_list = np.arange(np.sin(60*np.pi/180), np.sin((64 + 2/7)*np.pi/180), dx)[1:]
for theta in np.arcsin(x_list):
    xx, yy = n_fold(3, theta=theta)
    for n, (x, y) in enumerate(zip(xx, yy)):
        ax.plot(x, y, "C0", alpha=LOSS_BLUE**n)

ax.set_axis_off()
ax.set_aspect(1)
plt.show()
