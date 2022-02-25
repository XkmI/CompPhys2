from gpaw import GPAW , PW
import ase.io as io
from ase.build import bulk
from ase import Atoms
from ase.optimize import GPMin

babys_nth_calc = GPAW(xc ='PBE', mode=PW(450), kpts=(12, 12, 12), txt='calculation.txt')
all_that_glitters = bulk(name='Au')
the_sharpest_tool = bulk(name='Pt')
hey_now = bulk(name='Rh')

for a in [all_that_glitters, the_sharpest_tool, hey_now]:
    a.set_calculator(babys_nth_calc)
    dyn = GPMin(a, trajectory=f'{a}.traj', logfile=f'{a}.log')
    dyn.run(fmax=0.02, steps=100)
    a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()




#print(all_that_glitters)
#print(the_sharpest_tool)
#print(hey_now)


