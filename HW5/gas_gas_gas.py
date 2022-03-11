from ase import Atoms
from ase.build import molecule
from ase.optimize import GPMin, BFGS
from gpaw import GPAW, PW

o2 = molecule('O2', pbc=True)
co = molecule('CO', pbc=True)

o2.set_cell([12,12,12])
co.set_cell([12,12,12])
o2.center
co.center


babys_nth_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(1, 1, 1), spinpol=True, txt='calculation.txt')
babys_nplus1th_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(1, 1, 1), spinpol=False, txt='calculation.txt')

o2.set_calculator(babys_nth_calc)
relax = BFGS(o2, trajectory='o2.traj',logfile='o2.log')
relax.run(fmax=0.02)

co.set_calculator(babys_nplus1th_calc)
relax = BFGS(co, trajectory='co.traj',logfile='co.log')
relax.run(fmax=0.02)
