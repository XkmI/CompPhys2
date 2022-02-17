from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW, FermiDirac
import ase.io as io
from ase.optimize import GPMin

# The cluster we want to relax
atoms = io.read('christmas-tree.xyz')
# Define the calculator
calc = GPAW(nbands=10, h=0.25, txt='out.txt', occupations=FermiDirac(0.05), setups={'Na': '1'}, mode='lcao', basis='dzp')

atoms.set_calculator(calc)

dyn = GPMin(atoms, tratomsjectory='task8.traj', logfile='task8.log')
dyn.run(fmax=0.02, steps=100)
atoms.info['key_value_pairs']['raw_score'] = -atoms.get_potential_energy()
