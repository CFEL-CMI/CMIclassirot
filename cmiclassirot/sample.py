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


import numpy as np
import pyquaternion as quat
import scipy.constants
import pandas as pd



class Position(object):
    """Representation of a phase-space position of a molecule

    The current position can be directly accessed as `Position.angle`, which is a
    :class:`Quarternion` and `Position.velocity`, which is a 3-element `ndarray` -- but this will
    change in favor of them being implemented as property accessors in the future.

    .. todo:: Make `angle` and `velocity` a real property with accessors, so it get's implementation
    independent.

    .. todo:: Implement a 'raw' accessor for the internal data strucure, and epsecially a 'raw
    derivative'. The latter what should be set by Molecule._derivative and then be used through an
    accessor in the ODE solver... Need to figure out how to do this "efficiently".

    """

    def __init__(self, angle=quat.Quaternion(), velocity=np.zeros((3,)), t=0.):
        self.angle = angle
        self.velocity = velocity
        self.time = t



class Molecule(object):
    """Object to be manipulated

    The relevant parameters that describe the :class:`Molecule` properties are its principal moments
    of inertia and its polarizability tensor, both of which must be given in the principal axis of
    inertia system, i.e., :math:`a`, :math:`b`, :math:`c`.

    Furthermore, the molecule has a phase-space position, i.e., a :class:`Position object.

    """

    def __init__(self, *args, **kwargs):
        """Initialize a molecule from it's relevant paramters

        This can be called in two ways

        .. code-block:: python

            __init__(self, molecule, pos=None, T=None)
            __init__(self, I, P, angle=None, velocity=None, T=None, t=0.)

        :param molecule: Use first version of constructor: copy the specified molecule and possibly
            update its phase-space position

        :param pos: Determines position of new molecule
        * 'keep' or None specify to copy the phase-space position from the original Molecule
        * 'z' places the :class:`Molecule` along the laboratry :math:`z` axis with a angular
        velocity of 0
        * 'random' draws a random direction from a uniform angular distribution and an random
        angular velocity according to a 3D Maxwell distribution at temperature `temp`

        :param I: principal moments of inertia: these are the three diagonal elements of the in the
        inertial tensor in the principle axis of inertia frame (in SI units: kg * m**2)

        :param P: polarizabiity tensor in the inertial frame of the molecule (in SI units: ???)

        :param angle: initial position, as a Quaternion, if not provided a random direction id
        choosen

        :param velocity: initial angular frequency, if not provided a random vector drawn from three
        one-dimensional Maxell distribution for rotations about x, y, z is drawn

        :param T: Temperature of the thermal distribution from which to draw the initial momentum of
        the particle; typically this should not be done here but in a (to be provided) `Source`

        :param t: Time at which the molecule is created

        """
        if isinstance(args[0], Molecule):
            # first variant of constructor -- get args
            assert(len(args)==1)
            mol = args[0]

            # get keyword args
            pos = kwargs.get('pos', None)
            T = kwargs.get('T', 0)
            t = kwargs.get('t', 0)

            # construct object
            self._I = mol.I
            self._P = mol.P
            if pos == None or pos == 'keep':
                self.pos = mol.pos
            elif pos == 'z':
                self.pos = [Position()]
            elif pos == 'random':
                self.pos = [self._thermal_position(T, t)]
            elif isinstance(pos, Position(t=t)):
                self.pos = [pos]
            elif isinstance(pos, list):
                self.pos = pos
        elif len(args) == 0 or isinstance(args[0], np.ndarray):
            # second variant of constructor
            if len(args) >= 1:
                self._I = args[0]
            else:
                self._I = kwargs.get('I')
            if len(args) >= 2:
                self._P = args[1]
            else:
                self._P = kwargs.get('P')
            angle = kwargs.get('angle', None)
            velocity = kwargs.get('velocity', None)
            T = kwargs.get('T', None)
            t = kwargs.get('t', None)
            # set initila phase-space position
            self.pos = []
            pos = Position()
            if T is not None:
                pos = self._thermal_position(T)
            # overwrite thermal phase-space position by explicitely specified values
            if angle is not None:
                pos.angle = angle
            if velocity is not None:
                pos.velocity = velocity
            if t is not None:
                pass
                pos.time = t
            self.pos.append(pos)
        else:
            # something wrong!
            raise TypeError('Wrong argument types to Molecule.__init__')


    @property
    def I(self):
        return self._I


    @property
    def P(self):
        return self._P


    def acceleration(self, field, position, factor=1):
        """Calculate the molecule's angular acceleration at its current position in a field

        :param field: Field (amplitude, V/m) for which to calculate the acceleration


        .. todo:: Refactor code to get rid of the IFs.

        """
        q = position.angle
        v = position.velocity
        #I, P = self.rotate(q)  # uncomment this line incase you want to rotate
                                # the molecule not the field
        if self._I[2][2]<self._I[0][0]: # Particles initially oriented at the X-axis
            factor = 1.
        else: # Particles initially oriented along the Z-axis
            factor = -1.
        I, P = np.diagonal(self._I), self._P
        dipole = np.dot(P, field)
        torque = factor * np.cross(dipole, field)
        #a = np.dot(np.linalg.inv(I),   # if you decided to rotate the molecule
        #     (torque - np.cross(v, np.dot(I,v)))) #uncomment these lines and comment the next lines
        a = np.zeros(3)
        if I[0]!=0: # to deal with linear molecules
            a[0] = (torque[0] -(I[2] - I[1])*v[1]*v[2])/I[0]
        if I[1]!=0:
            a[1] = (torque[1] -(I[0] - I[2])*v[0]*v[2])/I[1]
        if I[2]!=0:
            a[2] = (torque[2] -(I[1] - I[0])*v[0]*v[1])/I[2]

        return a


    def _thermal_position(self, temp, t):
        """Calculate a random phase-space position for the specified temperature for generating
           random orientation of the particles, we rotate the molecules (initially on the Z-axis)
           around random axes on the x-y plane. The angle is calculated based on cosine destribution
           angle = acos(random.uniform(-1,1))

        .. todo:: this can be optimized, but let's write it down explicitely first

        """
        k = scipy.constants.Boltzmann
        velocity_0 = 0. # mean angular velocity is zero for every single dimension
        #velocity_sigma = np.sqrt(k*temp/np.diagonal(self._I)) # width of 1D Maxwell distribution
        #velocity = np.random.normal(velocity_0, velocity_sigma, 3)
        velocity_sigma = np.zeros(3)
        I = np.diagonal(self._I)
        if I[0]!=0:
            velocity_sigma[0] = np.sqrt(k*temp/I[0])
        if I[1]!=0:
            velocity_sigma[1] = np.sqrt(k*temp/I[1])
        if I[2]!=0:
            velocity_sigma[2] = np.sqrt(k*temp/I[2])

        velocity = np.random.normal(0., 1., 3)*velocity_sigma
        #print(velocity)
        #angle = acos(random.uniform(-1,1))
        #q = quat.Quaternion(axis=[random.uniform(-1,1),
        #    random.uniform(-1,1), 0], angle=angle)
        q = quat.Quaternion.random()
        return Position(q, velocity, t)


    def rotate(self, q):
        #R = q.inverse.rotation_matrix
        R = q.rotation_matrix.T
        I = np.dot(R, np.dot(self._I, R.T))
        P = np.dot(R, np.dot(self._P, R.T))
        return I, P



class Ensemble(object):
    """Source of randomly distributed Molecules

    This is also an iterator over :class:`Molecule`s

    """

    def __init__(self, size, molecule, T=0., t=0.):
        """Generate an ensemble  of |Molecule|s

        :param size: number of molecules in ensemble

        :param molecule:

        :param temperature:

        :param time: Time for the molecule to be create.

        .. todo:: add a constructor to create an ensemble from a data file

        """
        self.size = size
        self.molecules = []
        for i in range(self.size):
            self.molecules.append(Molecule(molecule, pos='random', T=T, t=t))
        self.index = self.size
        self.temperature = T
        self.time = t
        self.pulse = None


    @classmethod
    def FromFile(cls, name):
        cls._load(cls, name)
        return cls


    def __iter__(self):
        """Iterator method -- initialize"""
        return self


    def __next__(self):
        """Iterator method -- next"""
        if self.index == 0:
            raise StopIteration
        self.index -= 1
        return self.molecules[self.index]


    def _load(self, filename):
        self.molecules=[]
        with pd.HDFStore(filename) as store:
            tensors = store.get_storer('Molecule0').attrs.metadata
            self.size = len(store.keys())
            self.pulse = tensors['E']
            for i in store.keys():
               df = store[i]
               pos = []
               for j in df.values:
                   pos.append(Position(quat.Quaternion(j[:4]), j[4:7], t=j[7]))
               mol = Molecule(tensors['I'], tensors['P'])
               mol.pos = pos
               self.molecules.append(mol)
        self.index = self.size


    def _save(self, filename):
        df = pd.DataFrame(columns=['r', 'i', 'j', 'k',
                                   'omega_x', 'omega_y', 'omega_z', 'time'])
        store = pd.HDFStore(filename)
        for i in range(len(self.molecules)):
            df = pd.DataFrame(columns=df.columns)
            for k, j in enumerate(self.molecules[i].pos):
                df.loc[k] = list(np.concatenate((j.angle.elements, j.velocity)))+[j.time]
            store.put('Molecule'+str(i), df)
        metadata = {'I':self.molecules[0].I,'P':self.molecules[0].P,
                    'T':self.temperature, 'E': self.pulse}
        store.get_storer('Molecule0').attrs.metadata = metadata
        store.close()
