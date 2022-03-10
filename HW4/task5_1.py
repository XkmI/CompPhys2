from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum

lr = LrTDDFT('LrTDDFTresults.dat')
lr.diagonalize(energy_range=4)
photoabsorption_spectrum(lr,'spectrum_w.dat', width=0.06)