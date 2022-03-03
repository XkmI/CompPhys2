from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from gpaw import restart

_ , calc3 = restart('task1_backup.gpw')

calc3.get_eigenvalues()

dE = 6 # eV
lr = LrTDDFT(calc3, xc='LDA', restrict={'energy range': dE})
lr.diagonalize()
photoabsorption_spectrum(lr,'spectrum_w.dat',width=0.06)
lr.write('LrTDDFTresults.dat') # Save for task 2