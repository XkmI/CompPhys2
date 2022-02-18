from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW, FermiDirac
import ase.io as io
from ase.optimize import GPMin

# The cluster we want to relax
atoms = io.read('half-decahedron.xyz')
# Define the calculator
h0=0.25
#calc = GPAW(nbands=10, h=0.25, txt='out.txt', occupations=FermiDirac(0.05), setups={'Na': '1'}, mode='lcao', basis='dzp', xc='PBE')
calc = GPAW(mode='fd',h=h0,txt='out2.txt',symmetry={'point_group': False}, xc='PBE')
#calc = GPAW(mode='pw',h=h0,txt='out.txt',symmetry={'point_group': False})
atoms.set_calculator(calc)

dyn = GPMin(atoms, trajectory='half_fdpbe025.traj', logfile='half_fdpbe025.log')
dyn.run(fmax=0.02, steps=100)
atoms.info['key_value_pairs']['raw_score'] = -atoms.get_potential_energy()
