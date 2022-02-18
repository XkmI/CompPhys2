from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum

atoms = io.read('Na8.xyz')
atoms.center(vacuum=8.0) #add 8 ang vacuum

calc = GPAW(mode='fd', xc='LDA', setups={'Na': '1'}, h=0.3, nbands=0, txt='f.gpaw-out')

atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.write('task1_groundstate.gpw', mode='all') #saves the calculator, load it by calc=GPAW('babysFirstCalc(task1).gpw')

calc2 = GPAW('task1_groundstate.gpw')
calc2.set(nbands=110, convergence={'bands': -10}, fixdensity=True)

calc2.write('task1_backup.gpw', mode='all')
calc3 = GPAW('task1_backup.gpw')

calc3.get_eigenvalues()

dE = 6 # eV
lr = LrTDDFT(calc3, xc='LDA', restrict={'energy range': dE})
lr.diagonalize()
photoabsorption_spectrum(lr,'spectrum_w.dat',width=0.06)
lr.write('LrTDDFTresults.dat') # Save for task 2