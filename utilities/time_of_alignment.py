#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
import numpy as np
from cmiext import const

# physical constants
ε0 = const.vacuum_permittivity
c = const.speed_of_light # needed for unit conversion


# definition of an Au cylinder in its principal axis system
rho = 19300.    # density of gold (kg/m**3)
L = 2.e-9      # length of the cylinder (m)
R = 0.2e-9       # radius of the cylinder (m)
V = np.pi*L*R**2   # volume (m**3)
M = V*rho       # mass (kg)

# inertia tensor for a cylinder with the figure axis along z
# see ... (link to Wikipedia)
I = np.array([(0.083*M*L**2 + 0.25*M*R**2),
              (0.083*M*L**2 + 0.25*M*R**2),
               0.5*M*R**2])

### *** potential issue ***
###
### The following might not be the correct polarizability tensor -- which
### should be diagonal in the principal axis frame, not?
### Or is there something "special" about a *rod*?

# polarizability tensor for a cylinder  with the figure axis along z (m**3)
# obtained as polarizability volume from a ZENO calucation (angstrom**3),
# see ... (link to blogpost with input and documentation)
P = np.array([[ 50, 0.,  0.],
              [ 0., 50,  0.],
              [ 0., 0.,800.]])
# conversion to polariability in SI units (C m**2 / V)
# see https://en.wikipedia.org/wiki/Polarizability
P *= ε0 / 1e30


# definition of control field
# we assume a constant ("dc") field corresponding to a laser of intensity 10**10 W/cm**2
I0 = 1e14/2. # intensity (W/m**2)
# amplitude of the electric field (see https://www.rp-photonics.com/optical_intensity.html) in vacuum (n = 1)
A = np.sqrt((2*I0) / (c * ε0 * 1))
# applied electric field at 45º from Z in XZ plane
E = A * np.array((1,0,1))/np.sqrt(2)

# induced dipole: dot product between of polarizability and electric field
dipole = np.dot(P, E)
# torque: cross product between the dipole moment and E
torque = np.cross(dipole, E)
# angular acceleration (with center of mass at origin)
a = torque/I

# simple estimate of the time needed for the rod to rotate Δθ = 1/4 pi in space (Δθ = 1/2 a t**2)
t = np.sqrt(np.pi/4 / (np.abs(a)/2))

# print the result for the
print('Estimated time for 1/4 pi rotation about y (from θ=0..45º in the XZ plane): %.2g s' % t[1])
