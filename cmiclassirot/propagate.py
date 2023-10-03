# -*- coding: utf-8; fill-column: 100 -*-
#
# This file is part of CMIclassirot -- classical-physics rotational molecular-dynamics simulations
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE.md file
# for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy.integrate
import pyquaternion as quat
import multiprocessing as mp

from cmiclassirot.sample import Position


class Propagate(object):
    """Propagate an Ensemble in time"""

    def __init__(self, ensemble, field, timerange=(0,1e-9), dt_save=None):
        """Initialize propagator

        :param ensemble: :class:`Ensemble` with all |Molecule|s to be propagated

        :param field: :class:`Field` in which to propagate. This is a univariate function that
        returns the field strength (V/m) for any given time (s) in the :param:`timerange`.

        :param timerange: Time range over which to propagate, specified as tuple (t_init, t_final)

        :param dt_save: Timesteps for saving the current phase-space positions of the
        :param:`Molecule`s during propagation (default: 1 % of timerange period)

        """
        self.ensemble = ensemble
        self.field = field
        if dt_save:
            self.dt_save = dt_save
        else:
            self.dt_save = timerange[1] - timerange[0]
        self.t_range = timerange
        self.run()


    def run(self):
        """Propagate all |Molecule|s in the current |Field| over the current time range"""
        n_cpu = mp.cpu_count()
        print(f'Running on: {n_cpu} CPUs')
        self.ensemble.pulse = self.field
        pool = mp.Pool(n_cpu)
        results = pool.map(self._propagate, self.ensemble.molecules)
        pool.close()
        pool.join()
        for i in range(len(results)):
             self.ensemble.molecules[i] = results[i]
        return self.ensemble


    # @classmethod
    # def molecule(cls, mol, field, timerange, dt_save):
    #     cls.field = field
    #     cls.dt_save = dt_save
    #     cls.t_range = timerange
    #     cls._propagate(cls, mol)
    #     return mol


    def _propagate(self, molecule):
        """Propagate an individual molecule in the current field and over the current time range

        :param molecule: The :class:`Molecule` to propagate over the stored time-range and field; see `__init__`.

        .. note:: This method must be reentrant so we can run it for many molecules in parallel.

        .. todo:: Should use the modern approach, i.e., `scipy.integrate.solve_ivp`


        # create a copy store a field direction (which can the rotated)
        field = self.field
        derivative = np.empty((7,))
        rot = quat.Quaternion()
        # integrate over time
        scipy.integrate.solve_ivp(self._derivative,
                                  self.t_range,
                                  molecule.position,
                                  method='LSODA',
                                  t_eval=None,
                                  dense_output=False,
                                  events=None,
                                  vectorized=False,
                                  args=(derivative, molecule, field))
        """

        # initialize ode object
        integral = scipy.integrate.ode(self._derivative)
        # choose the integrator
        integral.set_integrator('dopri5', nsteps=10000)  # alternatively, use integral.set_integrator('lsoda')
        integral.set_initial_value(np.concatenate((molecule.pos[0].angle.elements,
                                                   molecule.pos[0].velocity)), self.t_range[0]).set_f_params(molecule)
        # integrate
        while integral.successful() and integral.t <= self.t_range[1]:
            integral.integrate(integral.t + self.dt_save)
            molecule.pos.append(Position(quat.Quaternion(integral.y[:4]), integral.y[4:7], t=integral.t))
        return molecule


    def _derivative(self, t, dpos, molecule):
        """Determine the derivative of a Molecule at a specific time.

        :param t: Time for which to calculate the derivative

        :param molecule: Molecule for which to calculate the derivative

        :param dpos: Derivative vector storage; if provided, this must be a 7-element ndarray which
        is used for storage of the derivative (which is returned in any case)

        :return: `ndarray` with the derivate of the current molecule position.

        .. note:: for performance reasons, i.e., avoiding the repeated creation of the derivative
        vector, this method expects the vector to be passed by the caller. This method then updates
        the elements and also returns the updated derivative vector.

        """
        # derivative of the posiiton is currently calculated under the assumption that pos.angle is a Quaternion

        q = quat.Quaternion(dpos[0:4])
        omega = dpos[4:7]
        # E = self.field(t) * self.field.rotate() # uncomment this if you want to rotate the molecule
        E = self.field(t) * self.field.rotate(q.inverse) # rotating the field
        position = Position(q, omega)
        dw = molecule.acceleration(E, position)
        dq = q.derivative(omega).elements
        return np.concatenate((dq,dw))
