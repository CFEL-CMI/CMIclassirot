#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of the CMI classical-rotation simulations

from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI classical-physics rotational molecular-dynamics simulations

This package provides facilities for the simulation of the rotation dynamics of molecules and particles in external
electric fields.

"""


setup(name="cmiclassirot",
      author              = "CFEL Controlled Molecule Imaging group",
      maintainer          = "CFEL Controlled Molecule Imaging group",
      maintainer_email    = "jochen.kuepper@cfel.de",
      url                 = "https://www.controlled-molecule-imaging.org/research/further_projects/software",
      description         = "CMI classical-physics rotational molecular-dynamics simulations",
      version             = "0.1.dev0",
      long_description    = long_description,
      license             = "GPL",
      packages            = ['cmiclassirot'],
      scripts             = ['scripts/cmiclassirot'],
      python_requires     = '>=3.5',
      install_requires    = ['numpy>=1.16.0',
                             'pyquaternion',
                             'scipy>=1.3.0'],
      )
