from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from gpaw import restart

atoms, calc2 = restart('task1_groundstate.gpw')
calc2.set(nbands=110, convergence={'bands': -10}, fixdensity=True)

calc2.write('task1_backup.gpw', mode='all')
