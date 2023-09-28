import numpy as np
from numpy import pi
from scipy.constants import h
 
# define a nanorod
rho = 19300.   # density *(kg/m**3)
L = 2.e-9      # length (m)
R = 2.e-10       # radius (m)
V = pi*L*R**2  # volume (m**3)
M = V*rho      # mass (kg)
 
# principal moments of inertia (see, e.g., Wikipedia)
I = np.array([(0.083*M*L**2 + 0.25*M*R**2),
              (0.083*M*L**2 + 0.25*M*R**2),
              0.5*M*R**2])
 
# calulate rotational constand B = \frac{h}{8 \pi^2 I} (Hz)
B = h / (8 * pi**2 * I)
 
# rotational periods (s)
rot_period = 1/(2*B)
 
print('rotational periods are calculated as ', rot_period, ' (s)')

