import numpy as np
from ase.units import Hartree
from gpaw.lrtddft import LrTDDFT
from helper.py import discrete_spectrum, fold, dump_data

calc=LrTDDFT('LrTDDFTresults.dat')
dump_data(calc, 'babysFirstDump.npz')
dump = np.load('babysFirstDump.npz')
K_pp = dump['K_pp']
dE_p = dump['ediff_p']

D_pp = np.diag(dE_p)
omega_pp = np.array([[0,D_pp],[-D_pp-2*K_pp,0]]) 