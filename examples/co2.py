from cmiclassirot.field import *
from cmiclassirot.sample import *
from scipy.constants import epsilon_0
from cmiclassirot.propagate import *
import numpy as np
P = 4 * np.pi * 1.e-30 * epsilon_0 * np.array([[ 1.93, 0.E2, 0.E2],
                                            [ 0.E2, 1.93, 0.E2],
                                            [ 0.E2, 0.E2, 4.03]])
#Polarizabilities are taken from 
#https://faculty.missouri.edu/~glaserr/vitpub/jp002927r.pdf

Iz = 0.
Ix = Iy = 7.14987936e-46 
I = np.array([[Ix ,0.,0.],
              [0., Iy,0.],
              [0., 0.,Iz]])
T = 298. 
sigma = 0.1e-12/ 2.355
# define a field
E = Field(amplitude=1.5e18, sigma=sigma, mean=0.25e-12)
n = 1500
# define a molecule
mol = Molecule(I, P)
# define an ensemble 
ensemble = Ensemble(n, mol, temperature=T)
# start simulation
p = Propagate(ensemble, E, timerange=(0.,1.5e-12), dt=1.e-14)
# save results
ensemble._save('CO2_150TWcm2_298K.h5')
