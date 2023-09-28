from cmiclassirot.field import *
from cmiclassirot.sample import *
from scipy.constants import epsilon_0
from cmiclassirot.propagate import *
import numpy as np

P = 1.e-30 * epsilon_0 * np.array([[ 3.06, 0.E2, 0.E2],
                                   [ 0.E2, 3.06, 0.E2],
                                   [ 0.E2, 0.E2, 7.46]])
Ix = Iy = 1.38493163e-45 
I = np.array([[Ix ,0.,0.],
              [0., Iy,0.],
              [0., 0.,0.]])
T = 0. 
sigma = 1.e-10 / 2.355
E = Field(amplitude=1.e14/2., sigma=sigma, mean=2.5e-10)
n = 50
mol = Molecule(I, P)
ensemble = Ensemble(n, mol, temperature=T)
p = Propagate(ensemble, E, timerange=(0.,1.e-9), dt=1.e-11)
ensemble._save('OCS_0.01TWcm2_0K.h5')
