# -*- coding: utf-8; fill-column: 100 -*-
#
# This file is part of CMIclassirot -- classical-physics rotational molecular-dynamics simulations
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file
# for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy import interpolate
from scipy.constants import c, epsilon_0
from pyquaternion import Quaternion

class Field(object):
    """Electric field class

    For the beginning, this only implements a simple Gaussian pulse (in field strengths!) centered
    at t=0 and along `Z`.

    .. todo:: Keep in mind that we want to be able to represent elliptically polaerized field by their envelope.

    """

    def __init__(self, amplitude=1.e15, sigma=1.e-9, mean=5.e-9, 
                 intensity=True, filename=None, dc=None):
        self.amplitude = amplitude
        self.sigma = sigma
        self.Ez = np.array([0., 0., 1.])
        self.from_file = None
        self.intensity = intensity
        self.mu = mean
        self.dc = dc
        if filename is not None:
            f = open(filename)
            E=[]
            t=[]
            for i in f:
                token = i.split()
                t.append(float(token[0]))
                E.append(float(token[1]))
            t = np.array(t)
            E = np.array(E)
            self.from_file = interpolate.interp1d(t,E)

    def __call__(self, t):
        """Field at time t

        :return: Field vector at time
        """
        if self.dc is not None:
            if t > self.dc[0] and t < self.dc[1]:
                return self.convert(self.dc[2])
            else:
                return 0

        if self.from_file is None:
            E = self.amplitude * np.exp(-0.5 * ((t-self.mu)/self.sigma)**2)
            if self.intensity:
                return self.convert(E)
            else:
                return E
        else:
            try:
                E = 1.e16*self.from_file(t)
                if self.intensity:
                    return self.convert(E)
                else:
                    return E
            except:
                return 0.
    
    def convert(self, I):
        """ Converts intensity to electric field"""
        return np.sqrt(2*I/(c * epsilon_0))

    def rotate(self, quaternion=Quaternion()):
        """Rotate field

        This is used in the actual propagation step, as it is easier/cheaper to rotate the field
        than to rotate the molecule

        :param rotation: Define the rotation to be applied to the field, as an :class:`Quarternion`

        """
        return quaternion.rotate(self.Ez)
