# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of CMIrotMD -- example calculation for a Au nanorod in a linear-plarization laser field
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
#
# This file shall eventually serve as a driver for the calculation and postprocessing
# for rotMD simulations


from cmiclassirot.field import *
from cmiclassirot.sample import *
from scipy.constants import epsilon_0
from cmiclassirot.propagate import *
import math
P = 1.e-30 * epsilon_0 * np.array([[ 5.646E5, 0.E2, 0.E2],
                                   [ 0.E2, 5.646E5, 0.E2],
                                   [ 0.E2, 0.E2, 4.860E6]])
nm = 1.e-9
rho    = 19.3e3 # density of gold (kg/m**3)
height = 10 * nm
radius = 2 * nm
volume = math.pi * radius**2 * height
mass   = volume*rho
Iz = 0.5*mass*radius**2
Ix = Iy = 1./12 * mass * (3*radius**2 + height**2)
I = np.array([[Ix ,0.,0.],
              [0., Iy,0.],
              [0., 0.,Iz]])

E = Field(amplitude=1.e15, sigma=1.e-9, mean=3.e-9, intensity=True)
n = 100
mol = Molecule(I, P)
ensemble = Ensemble(n, mol, temperature=20.)
p = Propagate(ensemble, E, timerange=(0.,15.e-9), dt=5.e-11)
ensemble._save('T298Knanorod.h5')
