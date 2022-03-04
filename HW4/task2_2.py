import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from helper import fold
'''
		Keywords for dumped quantities:
		'i_p':		Kohn-Sham state that gets excited from
		'a_p':		Kohn-Sham state that gets excited to
		'ediff_p':	Kohn-Sham eigenvalue differences for each KS excitation
		'fdiff_p':	Kohn_sham occupation differe for each KS excitation
		'mux_p':	x-component of dipole matrix elements
		'muy_p':	y-component of dipole matrix elements
		'muz_p':	z-component of dipole matrix elements
		'K_pp':		K-matrix
'''

dump = np.load('babysFirstDump.npz')

K_pp = dump['K_pp']
de_p = dump['ediff_p']

mu_p=np.zeros(3)
mu_p[0] = dump['mux_p']
mu_p[1] = dump['muy_p']
mu_p[2] = dump['muz_p']
df_p = dump['fdiff_p']
a_p = dump['a_p']
i_p = dump['i_p']

D_pp = np.diag(de_p)
omega_pp = np.block([[np.zeros(D_pp.shape),D_pp],[-D_pp-2*K_pp,np.zeros(D_pp.shape)]]) 
(X,F) = la.eig(omega_pp) #these are the discrete x vals

fi=0
for a in range(0,2):
    sum = 0
    for p in len(X):
        sum += mu_p[a][p] * np.sqrt(df_p[p]*de_p[p]) * F[p]
    fai = 2*np.abs(sum)*np.abs(sum)
    fi += fai
Y = fi #these are the discrete y vals

hartreeInEV = 27.211396641308
x = np.linspace(0,6/hartreeInEV)
y = fold(x,X,Y,0.06)

plt.plot(x,y)
plt.xlabel('Energy [Ha]')
plt.ylabel('Dipole moment [1/Ha]?')