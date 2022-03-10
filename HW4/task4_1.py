from ase import Atoms
from gpaw import GPAW
import ase.io as io
import gpaw.tddft as td

timestep=30      #as
iterations=1500  #1
power=1e-5       #???

xkick = [power, 0, 0]
ykick = [0, power, 0]
zkick = [0, 0, power]

td_calc_x = td.TDDFT('task1_groundstate.gpw')
td_calc_y = td.TDDFT('task1_groundstate.gpw')
td_calc_z = td.TDDFT('task1_groundstate.gpw')

#x-kick!
td_calc_x.absorption_kick(kick_strength=xkick)
td_calc_x.propagate(timestep, iterations, 'dipolemoment_x.dat', 'tddft_x.gpw')
td.photoabsorption_spectrum('dipolemoment_x.dat', 'x_spectrum.dat', width=0.06)

#y-kick!
td_calc_y.absorption_kick(kick_strength=ykick)
td_calc_y.propagate(timestep, iterations, 'dipolemoment_y.dat', 'tddft_y.gpw')
td.photoabsorption_spectrum('dipolemoment_y.dat', 'y_spectrum.dat', width=0.06)

#z-kick!
td_calc_z.absorption_kick(kick_strength=zkick)
td_calc_z.propagate(timestep, iterations, 'dipolemoment_z.dat', 'tddft_z.gpw')
td.photoabsorption_spectrum('dipolemoment_z.dat', 'z_spectrum.dat', width=0.06)

#No safety measures re:total time, as it allegedly is "sufficient".