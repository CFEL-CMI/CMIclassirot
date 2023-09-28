from cmiclassirot.sample import Ensemble
import matplotlib.pyplot as plt
import numpy as np
import sys

fname = sys.argv[1]
t=[]
cos2=[]
e = Ensemble.FromFile(fname)
r = np.zeros(len(e.molecules[0].pos))
E = e.pulse
print("Number of Molecules:", len(e.molecules))
for i in e.molecules:
    t=[]
    pulse = []
    for j in range(len(i.pos)):
     try:
      r[j] = r[j] + i.pos[j].angle.rotation_matrix[2][2]**2
      t.append(i.pos[j].time)
      pulse.append(E(i.pos[j].time))
     except:
      pass

r=r/len(e.molecules)
fig, ax = plt.subplots(2)
ax[0].plot(t, r)
ax[1].plot(t, pulse)
plt.xlabel("Time (s)")
ax[0].set_ylabel("Cos^2")
ax[1].set_ylabel("E (V/m)")
plt.tight_layout()
plt.show()
