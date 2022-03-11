from ase import Atoms
from ase.build import molecule
from ase.optimize import GPMin, BFGS
from gpaw import GPAW, PW
import ase.io as io

Pt_O = io.read('nPt_O.xyz')
#Pt_CO = io.read('nPt_CO.xyz')
Au_O = io.read('nAu_O.xyz')
#Au_CO = io.read('nAu_CO.xyz')
Rh_O = io.read('nRh_O.xyz')
#Rh_CO = io.read('nRh_CO.xyz')

calc1 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc2 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc3 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc4 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc5 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
calc6 = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')


Au_O.set_calculator(calc1)
relax = BFGS(Au_O, trajectory='nAu_O.traj',logfile='nAu_O.log')
relax.run(fmax=0.1)

#Au_CO.set_calculator(calc2)
#relax = BFGS(Au_CO, trajectory='nAu_CO.traj',logfile='nAu_CO.log')
#relax.run(fmax=0.1)

Pt_O.set_calculator(calc3)
relax = BFGS(Pt_O, trajectory='nPt_O.traj',logfile='nPt_O.log')
relax.run(fmax=0.1)

#Pt_CO.set_calculator(calc4)
#relax = BFGS(Pt_CO, trajectory='nPt_CO.traj',logfile='nPt_CO.log')
#relax.run(fmax=0.1)

Rh_O.set_calculator(calc5)
relax = BFGS(Rh_O, trajectory='nRh_O.traj',logfile='nRh_O.log')
relax.run(fmax=0.1)

#Rh_CO.set_calculator(calc6)
#relax = BFGS(Rh_CO, trajectory='nRh_CO.traj',logfile='nRh_CO.log')
#relax.run(fmax=0.1)
