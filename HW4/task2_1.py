import numpy as np
from gpaw.lrtddft import LrTDDFT
from helper.py import dump_data

calc=LrTDDFT('LrTDDFTresults.dat')
dump_data(calc, 'babysFirstDump.npz')
dump = np.load('babysFirstDump.npz')