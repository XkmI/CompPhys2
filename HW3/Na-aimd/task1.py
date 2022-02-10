#Vi behöver nog fler paket
import ase.io as io
import ase.md.npt as npt
from gpaw import GPAW
from ase.units import fs, kB 

atoms = io.read('cluster24_9176_Na.xyz')  

calc = GPAW(mode = 'lcao', xc = 'PBE', basis = 'dzp', symmetry = {'point_group': False}, charge = 1, txt = 'gpawOutput.gpaw-out')

atoms.set_calculator(calc)

dyn = npt.NPT(atoms, timestep = 0.5*fs, temperature = 350*kB, externalstress = 0, ttime = 20*fs, pfactor=None, logfile = 'babys_first_Nacluster24.log')

trajectory = io.Trajectory('babys_first_Nacluster24.traj', 'w', atoms)
dyn.attach(trajectory.write, interval=1)
dyn.run(4000)