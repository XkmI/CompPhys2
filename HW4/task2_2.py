import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

hartreeInEV = 27.211396641308

def fold(x_t, x_i, y_i, width):
    '''
        Convolutes each peak in the discrete spectrum
        with a Gaussian of chosen width.
        
        inputs:
        x_t:    vector with energies (i.e. linspace/np.arange)
        x_i:    stick spectrum energies
        y_i:    stick spectrum intensities
        width:  width of the Gaussian

        outputs:
        y_t:    convoluted spectrum, i.e. intensities
    '''
    def Gauss(x0):
        norm = 1.0 / (width * np.sqrt(2 * np.pi))
        return norm * np.exp(-0.5 * (x_t - x0)**2 / width**2)

    y_t = np.zeros_like(x_t)
    for x, y in zip(x_i, y_i):
        y_t += y * Gauss(x)
    return y_t

'''
		Keywords for dumped quantities:
		'i_p':		Kohn-Sham state that gets excited from
		'a_p':		Kohn-Sham state that gets excited to
		'ediff_p':	Kohn-Sham eigenvalue differences for each KS excitation
		'fdiff_p':	Kohn_sham occupation differences for each KS excitation
		'mux_p':	x-component of dipole matrix elements
		'muy_p':	y-component of dipole matrix elements
		'muz_p':	z-component of dipole matrix elements
		'K_pp':		K-matrix
'''

dump = np.load('babysFirstDump.npz')

K_pp = dump['K_pp']
de_p = dump['ediff_p']

#mu_p=np.zeros(3)
mux_p = dump['mux_p']
muy_p = dump['muy_p']
muz_p = dump['muz_p']
df_p = dump['fdiff_p']
a_p = dump['a_p']
i_p = dump['i_p']

D_pp = np.diag(de_p)
omega_pp = np.matmul(np.matmul(np.sqrt(D_pp),(D_pp)),np.sqrt(D_pp)) # +2*K_pp
(eigs,F) = la.eig(omega_pp) 
X = np.sqrt(eigs) #these are the discrete x vals
#print(omega_pp)
mu = np.array([mux_p, muy_p, muz_p])
fI = np.zeros(mu.shape)
#print(F[:,10:11])
for a in range(mu.shape[0]):
	for i in range(mu.shape[1]):
		sum = 0
		for p in range(mu.shape[1]):
			sum += mu[a][p] * np.sqrt(df_p[p]*de_p[p]) * F[p,i]
		fI[a][i] = 2 * np.abs(sum)**2

Y = np.sum(fI,0)/3 #these are the discrete y vals
x = np.linspace(0/hartreeInEV,6.4/hartreeInEV,500) # 6/hartreeInEV
#plt.plot(x,Y)
y = fold(x,X,Y,0.06/hartreeInEV)

plt.plot(x*hartreeInEV,y/hartreeInEV)
plt.xlabel('Energy [eV]')
plt.ylabel('Dipole moment [1/eV]')
plt.show()