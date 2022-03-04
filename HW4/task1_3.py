from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from gpaw import restart

calc3 = GPAW('task1_groundstate110.gpw')
calc3.get_eigenvalues()

dE = 6 # eV
lr = LrTDDFT(calc3, xc='LDA', energy_range=dE)
lr.diagonalize()
photoabsorption_spectrum(lr,'spectrum_w.dat', width=0.06)
lr.write('LrTDDFTresults.dat') # Save for task 2