# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of the CMIclassirot classical-rotation simulations

"""Example CMIclassirot calculation for the impulsive alignment of cold OCS

Laser alignment of OCS sample at 2 K with a FWHM = 500 fs Gaussian pulse with an peak intensity of $10^{13}$ W/cm$^2$

@author: Jochen Küpper <jochen.kuepper@cfel.de>

"""

from math import log, sqrt
import numpy as np
import scipy
from scipy.constants import c, epsilon_0, h, pi, physical_constants
a0 = physical_constants['Bohr radius'][0]

from cmiclassirot.field import *
from cmiclassirot.sample import *


# define timerange of integration and stepsize for saving the distribution
timerange = (-10e-12, 40e-12)   # simulate from -10...40 ps
dt_save   = 250e-15              # and save positions every 250 fs



# define the alignment field
# We use a Gaussian temporal profile with a FWHM = 500 fs, centered at 0,
# and a peak intensity of 10^13 W/cm**2
field = Field(peak_intensity=1e13 * 1e4, FWHM=500e-15, t_peak=0.)



# define the polarizability tensor of OCS in atomic units (a_0^3) [Holmegaard et al., Nature Physics 6, 428 (2010)]
P = np.array([[26.15, 0.,    0.],
              [0.,    26.15, 0.],
              [0.,    0.,    50.72]])
# and convert from polarizability volume (a_0^3) to polarizability (C m^2 V^-1)
P *= 4 * pi * epsilon_0 * a0**3

# define the inertial tensor of OCS (kg m^2)
#
# OCS rotational constant is B = 0.20286 cm-1 [Herzberg, Electronic spectra and electronic structure of polyatomic
# molecules, (van Nostrand, New York, 1966)]
B = 0.20286 * 1e-2
I = h / (8 * pi**2 * c) * B
I = np.array([[I, 0, 0],
              [0, I, 0],
              [0, 0, 0]])

# define Molecule based on OCS values
mol = Molecule(I, P, t=timerange[0])

# Perform calculation for initial temperature of 2 K
T = 2.

# we calculate the alignment for an ensemble of 10,000 randomly sampled molecules
n = 10000
ensemble = Ensemble(n, mol, T=T, t=timerange[0])
