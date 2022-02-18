from ase import Atoms
from gpaw import GPAW
import ase.io as io

atoms = io.read('Na8.xyz')
atoms.center(vacuum=8.0) #add 8 ang vacuum

calc = GPAW(mode='fd', xc='LDA', setups={'Na': '1'}, h=0.3, nbands=0, txt='f.gpaw-out')

atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.write('babysFirstCalc(task1).gpw', 'all') #saves the calculator, load it by calc=GPAW('babysFirstCalc(task1).gpw')

calc.set(nbands=110, convergence={'bands': -10}, fixdensity=True)