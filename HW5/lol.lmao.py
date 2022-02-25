from gpaw import GPAW , PW
import ase.io as io
from ase.build import bulk
from ase import Atoms
from ase.optimize import GPMin

babys_nth_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(12, 12, 12), txt='calculation.txt')
babys_nplus1th_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(12, 12, 12), txt='calculation.txt')
babys_nplus2th_calc = GPAW(xc ='PBE', mode=PW(450), kpts =(12, 12, 12), txt='calculation.txt')
all_that_glitters = bulk(name='Au')
the_sharpest_tool = bulk(name='Pt')
hey_now = bulk(name='Rh')


all_that_glitters.set_calculator(babys_nth_calc)
dyn = GPMin(all_that_glitters, trajectory='all_that_glitters.traj', logfile='all_that_glitters.log')
dyn.run(fmax=0.02, steps=100)
all_that_glitters.info['raw_score'] = -all_that_glitters.get_potential_energy()

the_sharpest_tool.set_calculator(babys_nplus1th_calc)
dyn = GPMin(the_sharpest_tool, trajectory='the_sharpest_tool.traj', logfile='the_sharpest_tool.log')
dyn.run(fmax=0.02, steps=100)
the_sharpest_tool.info['raw_score'] = -the_sharpest_tool.get_potential_energy()

hey_now.set_calculator(babys_nplus2th_calc)
dyn = GPMin(hey_now, trajectory='hey_now.traj', logfile='hey_now.log')
dyn.run(fmax=0.02, steps=100)
hey_now.info['raw_score'] = -hey_now.get_potential_energy()



#print(all_that_glitters)
#print(the_sharpest_tool)
#print(hey_now)


