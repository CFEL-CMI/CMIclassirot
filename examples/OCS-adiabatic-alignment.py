# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of the CMIclassirot classical-rotation simulations

"""Example CMIclassirot calculation: long-pulse alignment of cold OCS with a FWHM = 1 ns Gaussian pulse

This example uses the rotational constanf of OCS from [Herzberg, Electronic spectra and electronic structure of
polyatomic molecules, (van Nostrand, New York, 1966)] and the polarizability tensor from [Holmegaard et al., Nature
Physics 6, 428 (2010)]

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
timerange = (-2.5e-9, 2.5e-9)  # simulate from -2.5...2.5 ns
dt_save   = 50e-12             # and save positions every 50 ps


# define the alignment field
# We use a Gaussian temporal profile with a FWHM = 1 ns
# and a peak intensity of 10**12 W/cm**2
field = Field(peak_intensity=1e12 * 1e4, FWHM=1e-9, t_peak=0)


# define the polarizability tensor of OCS in atomic units (a_0^3) [Holmegaard et al., Nature Physics 6, 428 (2010)]
P = np.array([[26.15, 0.,    0.],
              [0.,    26.15, 0.],
              [0.,    0.,    50.72]])
# and convert from polarizability volume (a_0^3) to polarizability (C m^2 V^-1)
P *= 4 * pi * epsilon_0 * a0**3

# define the inertial tensor of OCS (kg m^2)
#
# rotational constant is B = 0.20286 cm-1 [Herzberg, Electronic spectra and electronic structure of polyatomic
# molecules, (van Nostrand, New York, 1966)]
B = 0.20286 * 1e-2
I = h / (8 * pi**2 * c) * B
I = np.array([[I, 0, 0],
              [0, I, 0],
              [0, 0, 0]])

# define Molecule based on OCS values
mol = Molecule(I, P, t=timerange[0])

# Perform calculation for initial temperature of 0 K
T = 0.

# we calculate the alignment for an ensemble of 1000 randomly sampled molecules
n = 1000
ensemble = Ensemble(n, mol, T=T, t=timerange[0])
