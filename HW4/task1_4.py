import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('spectrum_w.dat',skiprows=4)
y=data[:,1]
x=data[:,0]
plt.plot(x,y)
plt.xlabel('Energy [eV]')
plt.ylabel('Dipole moment [1/eV]')
plt.show()