# -*- coding: utf-8; fill-column: 100 -*-
#
# This file is part of CMIclassirot -- classical-physics rotational molecular-dynamics simulations
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# If you use this programm for scientific work, you must correctly reference it; see the LICENSE.md
# file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

__doc__ = """Running this file gives mane warnings -- fix.
.. todo:: SyntaxWarning: "is" with a literal. Did you mean "=="?
"""

from cmiclassirot.sample import *
from cmiclassirot.field import *
from math import pi
import unittest
from pyquaternion import Quaternion
from scipy.constants import c, epsilon_0

# moments of inertia of solid cylinder of gold
# see https://en.wikipedia.org/wiki/List_of_moments_of_inertia
rho    = 19.3e3 # density of gold (kg/m**3)
height = 20e-9 # m
radius = 1e-9 # m
volume = pi * radius**2 * height
mass   = volume*rho

# inertia tensor in the PA frame
I = np.array([1./12 * mass * (3*radius**2 + height**2),
              1./12 * mass * (3*radius**2 + height**2),
              0.5*mass*radius**2])

# polarizability tensor of a perfectly conducting nanorod of 20 nm height and 2 nm diameter
P = np.array([[  1.1e-26, 0.,       0],
              [  0.00e0,  1.1e-26,  0],
              [  0,       0,  2.5e-25]])

T = 0.4 #K
# testing the first version of the constructor of the molecule class
class TestCMIclassirot(unittest.TestCase):

    def test_molecule_constructor_first_version(self):
        """
        Testing the first version of the constructor of the
        molecule class.
        """

        m1 = Molecule(I, P)
        m2 = Molecule(m1)
        self.assertEqual(m1.I.tolist(), m2.I.tolist())
        self.assertEqual(m1.P.tolist(), m2.P.tolist())
        m2 = Molecule(m1, pos='random')
        self.assertIsNotNone(m2.pos[0].angle)
        self.assertIsNotNone(m2.pos[0].velocity)
        velocity = np.zeros((3,))
        angle = Quaternion()
        m2 = Molecule(m1, pos='z')
        self.assertEqual(m2.pos[0].velocity.tolist(), velocity.tolist())
        self.assertEqual(m2.pos[0].angle.angle, angle.angle)
        pos = Position(angle, velocity)
        m2 = Molecule(m1, pos=pos)
        self.assertEqual(m2.pos[0].velocity.tolist(), velocity.tolist())
        self.assertEqual(m2.pos[0].angle.angle, angle.angle)

    def test_molecule_constructor_second_version(self):
        """
        Testing the second version of the constructor of the
        molecule class.
        """

        m = Molecule(I, P, temp=T)
        self.assertEqual(m.I.tolist(), I.tolist())
        self.assertEqual(m.P.tolist(), P.tolist())
        self.assertIsNotNone(m.pos[0].angle)
        self.assertIsNotNone(m.pos[0].velocity)

    def test_assigning_speed_angles(self):
        """
        test passing values for angle and velocities
        """

        velocity = np.array([3e6, 2e6, -1.e7])
        angle = Quaternion()
        m = Molecule(I, P, temp=T, velocity=velocity, angle=angle)
        self.assertEqual(m.pos[0].velocity.tolist(), velocity.tolist())
        self.assertEqual(m.pos[0].angle.angle, angle.angle)

    def test_Ensemble(self):
       """
       test generating ensemble of molecules
       """
       n = 100
       mol = Molecule(I, P)
       ensemble = Ensemble(n, mol)
       self.assertEqual(ensemble.size, n)
       self.assertEqual(ensemble.size, ensemble.index)

    def test_save_load_ensemble(self):
        """ 
        testing save and load of ensemble 
        to and from file"""
        n = 100
        mol = Molecule(I, P)
        ensemble = Ensemble(n, mol)
        ensemble._save('Data.h5')
        ensemble_from_file = Ensemble.FromFile('Data.h5')
        self.assertEqual(len(ensemble.molecules), len(ensemble_from_file.molecules))

    def test_intensity_to_field(self):
        """testing the conversion from intenisty 
           to electric field. if the intensity = c*epsilon_0
           then the field should equal 1 V/d"""

        A = c * epsilon_0/2.
        E = Field(amplitude=A, sigma=5e-9, mean=0., intensity=True)
        self.assertEqual(E(0.)[2], 1.)

if __name__ == '__main__':
    unittest.main()
