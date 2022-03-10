import numpy as np
import matplotlib.pyplot as plt

xdata = np.loadtxt('x_spectrum.dat',skiprows=4)
yx=xdata[:,1]
xx=xdata[:,0]
ydata = np.loadtxt('y_spectrum.dat',skiprows=4)
yy=ydata[:,1]
xy=ydata[:,0]
zdata = np.loadtxt('z_spectrum.dat',skiprows=4)
yz=zdata[:,1]
xz=zdata[:,0]

plt.plot(xx,yx,label="Field $\\propto \\hat{x}$ ")
plt.plot(xy,yy,label="Field $\\propto \\hat{y}$ ")
#plt.plot(xz,yz,label="Field $\\propto \\hat{z}$ ")

plt.xlabel('Energy [eV]')
plt.ylabel('Dipole moment [1/eV]')
plt.legend()
plt.axis([0, 7, -3, 30])
plt.show()