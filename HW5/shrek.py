import numpy as np
from ase import Atoms
from ase.build import bulk
# from ase.io import write
# from ase.io.trajectory import PickleTrajectory
from ase.io.trajectory import Trajectory
# from ase.optimize.bfgs import BFGS
from gpaw import GPAW
from gpaw import PW
# from gpaw import FermiDirac

atomNames = ['Au', 'Pt', 'Rh']
latPms = [4.15, 3.95, 3.85]

for i in range(3):
    all_that_glitters = bulk(atomNames[i], 'fcc', a=latPms[i])

    cell = all_that_glitters.get_cell()
    #traj = PickleTrajectory(atomNames[i] + '.traj', 'w')

    calc=GPAW(mode=PW(450),                # Energycutoff for planewaves [eV]
              h=0.2,                       # The distance between gridpoints AA^-1
              xc='PBE',                    # xc-functional
              kpts=(12,12,12),             # number of k-points
              # occupations=FermiDirac(0.1), # Fermi temperature [eV]
              txt=atomNames[i]+'_bulk_a.txt')
            
    all_that_glitters.set_calculator(calc)

    for x in np.linspace(0.95, 1.05, 50):
        all_that_glitters.set_cell(cell * x, scale_atoms=True)
        all_that_glitters.info['raw_score'] = all_that_glitters.get_potential_energy()
        traj = Trajectory(atomNames[i] + '.traj', mode='a', atoms=all_that_glitters)
        traj.write()