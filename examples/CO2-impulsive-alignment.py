# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of the CMIclassirot classical-rotation simulations

"""Example CMIclassirot calculation for the impulsive alignment of warm CO2

Laser alignment of CO2 sample at room temperature 296 K with a FWHM = 100 fs Gaussian pulse with an peak intensity of
:math:`$1.5\cdot10^{14}$ W/cm$^2$`

Compare results to Fig. 1 of [[J. M. Hartmann and C. Boulet, J. Chem. Phys. 136, 184302
(2012)]](https://dx.doi.org/10.1063/1.4705264).

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
timerange = (-1e-12, 3e-12)      # simulate from -1...3 ps
dt_save   = 10e-15               # and save positions every 50 fs



# define the alignment field
# We use a Gaussian temporal profile with a FWHM = 100 fs, centered at 0,
# and a peak intensity of 1.5 * 10^14 W/cm**2
field = Field(peak_intensity=1.5e14 * 1e4, FWHM=100e-15, t_peak=0.)



# define the polarizability tensor of CO2 in atomic units (a_0^3) [Bridge, N. J.; Buckingham, A. D. Proc. R. Soc. A,
# 295, 33 (1966)]
P = np.array([[13.018,  0.,     0.],
              [ 0.,    13.018,  0.],
              [ 0.,     0.,    27.250]])
# and convert from polarizability volume (a_0^3) to polarizability (C m^2 V^-1)
P *= 4 * pi * epsilon_0 * a0**3

# define the inertial tensor of OCS (kg m^2)
#
# OCS rotational constant is B = 0.20286 cm-1 [Herzberg, Electronic spectra and electronic structure of polyatomic
# molecules, (van Nostrand, New York, 1966)]
B = 0.39021 * 1e-2
I = h / (8 * pi**2 * c) * B
I = np.array([[I, 0, 0],
              [0, I, 0],
              [0, 0, 0]])

# define Molecule based on OCS values
mol = Molecule(I, P, t=timerange[0])

# Perform calculation for initial temperature of 298 K (room temeprature)
T = 296.

# we calculate the alignment for an ensemble of 10,000 randomly sampled molecules
n = 10000
ensemble = Ensemble(n, mol, T=T, t=timerange[0])
