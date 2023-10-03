#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of the CMI classical-rotation simulations

from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI classical-physics rotational molecular-dynamics simulations (CMIclassirot)

This package provides facilities for the simulation of the rotation dynamics of molecules and particles in external
electric fields.

"""

name = 'CMIclassirot'
release = '0.9.1'
version = '0.9.1'

setup(name=name,
    author              = "CFEL Controlled Molecule Imaging group",
    maintainer          = "CFEL Controlled Molecule Imaging group",
    maintainer_email    = "jochen.kuepper@cfel.de",
    url                 = "https://www.controlled-molecule-imaging.org/research/further_projects/software",
    license             = "GPLv3",
    description         = "CMI classical-physics rotational molecular-dynamics simulations",
    long_description    = long_description,
    version             = version,
    packages            = ['cmiclassirot'],
    scripts             = ['bin/cmiclassirot',
                           'bin/cmiclassirot-plot'],
    python_requires     = '>=3.9',
    install_requires    = ['numpy>=1.16.0',
                           'pyquaternion',
                           'scipy>=1.3.0',
                           'tables'],
    )
