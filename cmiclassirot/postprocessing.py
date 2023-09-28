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


class Expectation(object):
    """Calculate expectation values of a probability density function given as an :class:`Ensemble`

    .. todo:: implement

    """

    def __init__(self, Ensemble):
        pass

    def __call__(self):
        pass


class cos2theta(Expectation):
    pass


class cos2theta_2D(Expectation):
    pass



class ProbabilityGraphics(object):
    """Create a graphical representation of a probability density function

    This Object requires a (variable dimension) probability density function and provides various
    graphical outputs, such as actual 3D images/VRML graphs, contour plots of projections, etc.

    .. todo:: implement

    """

    def __init__(self):
        pass



class Projection(object):
    """Project a probability distribution onto a plane or line

    This object takes a 3D PDF and creates a lower-dimension PDF by projection onto an arbitrary
    plane or line.

    .. todo:: Discuss how this blends into :class:`Ensemble` and the funtionality in this file: For
    instance, will we just create a new Ensemble with certain restrictions or in fact create an
    independent data structure?

    .. todo:: implement

    """

    def __init__(self):
        pass
