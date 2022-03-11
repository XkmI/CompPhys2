from ase import Atoms
from ase.build import fcc111
from ase.optimize import GPMin
from gpaw import GPAW, PW

Rhodium_ingot = fcc111('Rh', size=(3,3,3), vacuum=6.0, a=3.84)
Gold_ingot = fcc111('Au', size=(3,3,3), vacuum=6.0, a=4.18)
Platinum_ingot = fcc111('Pt', size=(3,3,3), vacuum=6.0, a=3.97)

babys_nth_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
babys_nplus1th_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')
babys_nplus2th_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(4, 4, 1), txt='calculation.txt')

Rhodium_ingot.set_calculator(babys_nth_calc)
dyn = GPMin(Rhodium_ingot, trajectory='newRhodium_ingot.traj', logfile='newRhodium_ingot.log')
dyn.run(fmax=0.02, steps=100)
Rhodium_ingot.info['raw_score'] = -Rhodium_ingot.get_potential_energy()

Gold_ingot.set_calculator(babys_nplus1th_calc)
dyn = GPMin(Gold_ingot, trajectory='newGold_ingot.traj', logfile='newGold_ingot.log')
dyn.run(fmax=0.02, steps=100)
Gold_ingot.info['raw_score'] = -Gold_ingot.get_potential_energy()

Platinum_ingot.set_calculator(babys_nplus2th_calc)
dyn = GPMin(Platinum_ingot, trajectory='newPlatinum_ingot.traj', logfile='newPlatinum_ingot.log')
dyn.run(fmax=0.02, steps=100)
Platinum_ingot.info['raw_score'] = -Platinum_ingot.get_potential_energy()