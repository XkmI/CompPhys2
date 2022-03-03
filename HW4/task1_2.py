from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
from gpaw import restart

atoms = io.read('Na8.xyz')
atoms.center(vacuum=8.0) #add 8 ang vacuum

calc2 = GPAW('task1_groundstate.gpw')
calc2.set(nbands=110, convergence={'bands': -10}, fixdensity=True)

atoms.set_calculator(calc2)
atoms.get_potential_energy()

calc2.write('task1_groundstate110.gpw', mode='all')
