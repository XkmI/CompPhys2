from ase import Atoms
from ase.build import molecule
from ase.optimize import GPMin, BFGS
from gpaw import GPAW, PW
import ase.io as io

#Pt_O = io.read('Pt_O.xyz')
Pt_CO = io.read('Pt_CO.xyz')
#Au_O = io.read('Au_O.xyz')
Au_CO = io.read('Au_CO.xyz')
#Rh_O = io.read('Rh_O.xyz')
Rh_CO = io.read('Rh_CO.xyz')

#calc1 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc2 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
#calc3 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc4 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
#calc5 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc6 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')


#Au_O.set_calculator(calc1)
#relax = BFGS(Au_O, trajectory='Au_O.traj',logfile='Au_O.log')
#relax.run(fmax=0.1)

Au_CO.set_calculator(calc2)
relax = BFGS(Au_CO, trajectory='Au_CO.traj',logfile='Au_CO.log')
relax.run(fmax=0.1)

#Pt_O.set_calculator(calc3)
#relax = BFGS(Pt_O, trajectory='Pt_O.traj',logfile='Pt_O.log')
#relax.run(fmax=0.1)

Pt_CO.set_calculator(calc4)
relax = BFGS(Pt_CO, trajectory='Pt_CO.traj',logfile='Pt_CO.log')
relax.run(fmax=0.1)

#Rh_O.set_calculator(calc5)
#relax = BFGS(Rh_O, trajectory='Rh_O.traj',logfile='Rh_O.log')
#relax.run(fmax=0.1)

Rh_CO.set_calculator(calc6)
relax = BFGS(Rh_CO, trajectory='Rh_CO.traj',logfile='Rh_CO.log')
relax.run(fmax=0.1)
