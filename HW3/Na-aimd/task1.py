import time
start_time = time.time()  

#Vi behöver nog fler paket
from ase.units import fs, kB 
from gpaw import GPAW

atoms = read('cluster24.log')  #känns bootleg

calc = GPAW(
    mode = 'lcao',
    xc = 'PBE',
    basis = 'dzp'
    symmetry = {'point_group': False},
    charge = 1,
    txt = 'gpawOutput.gpaw-out',
)

atoms.set_calculator(calc)

dyn = NPT(
    temperature = 350*kB,
    timestep = 0.5*fs,
    ttime = 20*fs,
    logfile = 'mdOutput.log',
)

trajectory = Trajectory('someDynamics.traj', 'w', atoms)  #känns bootleg
dyn.attach(trajectory.write, interval=1)
dyn.run(10)