from ase import Atoms
from gpaw import GPAW
import ase.io as io
from gpaw.tddft import *

timestep=30      #as
iterations=1500  #1
power=1e-5       #???

xkick = [power, 0, 0]
ykick = [0, power, 0]
zkick = [0, 0, power]

td_calc_x = TDDFT('task1_groundstate.gpw')
td_calc_y = TDDFT('task1_groundstate.gpw')
td_calc_z = TDDFT('task1_groundstate.gpw')

#x-kick!
DipoleMomentWriter(td_calc_x, 'dipolemoment_x.dat')
RestartFileWriter(td_calc_x, 'tddft_x.gpw')
td_calc_x.absorption_kick(kick_strength=xkick)
td_calc_x.propagate(timestep, iterations)
td_calc_x.write('tddft_x.gpw', mode='all')
photoabsorption_spectrum('dipolemoment_x.dat', 'x_spectrum.dat', width=0.06)

#y-kick!
DipoleMomentWriter(td_calc_y, 'dipolemoment_y.dat')
RestartFileWriter(td_calc_y, 'tddft_y.gpw')
td_calc_y.absorption_kick(kick_strength=ykick)
td_calc_y.propagate(timestep, iterations)
td_calc_y.write('tddft_y.gpw', mode='all')
photoabsorption_spectrum('dipolemoment_y.dat', 'y_spectrum.dat', width=0.06)

#z-kick!
DipoleMomentWriter(td_calc_z, 'dipolemoment_z.dat')
RestartFileWriter(td_calc_z, 'tddft_z.gpw')
td_calc_z.absorption_kick(kick_strength=zkick)
td_calc_z.propagate(timestep, iterations)
td_calc_z.write('tddft_z.gpw', mode='all')
photoabsorption_spectrum('dipolemoment_z.dat', 'z_spectrum.dat', width=0.06)

#No safety measures re:total time, as it allegedly is "sufficient".