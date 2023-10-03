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

    For the beginning, this only implements a simple Gaussian pulse (in amplitude) centered at t=0
    and along `Z`.

    Must specify one of `peak_amplitude`, `peak_intensity`, or `fieldname` to provide an
    ac-electric-field envelope. In the formre two cases, one must also specify the width as `FWHM`
    in amplitude or intensity, respectively, and the time of the `peak` of the laser pulse. In the
    latter case, `file_content` needs to be specified.

    :param peak_amplitude: Peak amplitude of a synthetic Gaussian laser pulse (V/m)

    :param peak_intensity: Peak intensity of a synthetic Gaussian laser pulse (W/cm^2)

    :param t_peak: Time of the laser-pulse peak (s)

    :param FWHM: Full-width-at-half-maximum width of laser pulse in same units as filed
    specification, i.e., :param`peak_amplitude` or :param:`peak_intensity` (s)

    :param filename: Name of file to read numerically-specified field

    :param file_content: Specify if file contains `amplitude` (V/m) or `intensity` (W/c^2) values.
    This functionality is currently not implemented and the file is always expwected to contain
    amplitudes. Please improve.

    .. todo:: Keep in mind that we want to be able to represent elliptically polaerized field by their envelope.

    .. todo:: document class, constructor, and methods

    .. todo:: Implement mixed field handling, here this is always providing (DC, AC) tuples on
    __call__; maybe should be done via a derived class MixedField for performance reasons.

    """
    def __init__(self, peak_amplitude=None, peak_intensity=None, t_peak=0., FWHM=10.e-9,
                 filename=None, file_content='intensity'):
        # store field internally as an field apmplitude
        assert (peak_amplitude or peak_intensity or filename) # at least one is not None
        if peak_intensity:
            self.peak_amplitude = Field.intensity2amplitude(peak_intensity)
            self.sigma = FWHM / 2*np.sqrt(2*np.log(2))
            self.t_peak = t_peak
        elif peak_amplitude:
            self.peak_amplitude = peak_amplitude
            self.sigma = FWHM / 2*np.sqrt(2*np.log(2))
            self.t_peak = t_peak
        else:
            self.peak_amplitude = None
            f = open(filename)
            E=[]
            t=[]
            for i in f:
                token = i.split()
                t.append(float(token[0]))
                E.append(float(token[1]))
            t = np.array(t)
            E = np.array(E)
            self.amplitude = interpolate.interp1d(t,E)
        # initially the field is alog :math:`Z`
        self.Ez = np.array([0., 0., 1.])


    def __call__(self, t):
        """Field at time t

        :return: Field amplitude vector at time :math:`t`

        """
        if self.peak_amplitude:
            return self.peak_amplitude * np.exp(-0.5 * ((t - self.t_peak)/self.sigma)**2)
        else:
            return self.amplitude(t)


    @staticmethod
    def intensity2amplitude(I):
        """Converts intensity (in W/m**2) to electric field (in V/m)"""
        return np.sqrt(2 * I / (c * epsilon_0))


    @staticmethod
    def amplitude2intensity(A):
        """Converts electric field amplitude (in V/m) to intensity (in W/m**2)"""
        return  c * epsilon_0 / 2 * A**2


    def rotate(self, quaternion=Quaternion()):
        """Rotate field

        This is used in the actual propagation step, as it is easier/cheaper to rotate the field
        than to rotate the molecule

        :param rotation: Define the rotation to be applied to the field, as an :class:`Quarternion`

        """
        return quaternion.rotate(self.Ez)
