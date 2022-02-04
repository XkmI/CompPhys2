from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt


d=1.1
#molekyl = atoms('n st atom X', [xyz],[xyz].....)
molecule6 = Atoms('6Na', [0,0,0],[0,0,d],[0,d,0],[d,0,0],[d,d,0],[0,d,d])
molecule7 = Atoms('7Na', [0,0,0],[0,0,d],[0,d,0],[d,0,0],[d,d,0],[0,d,d],[d,0,d])
molecule8 = Atoms('8Na', [0,0,0],[0,0,d],[0,d,0],[d,0,0],[d,d,0],[0,d,d],[d,0,d],[d,d,d]) #Simple noobic
molecule6.calc = EMT()
molecule7.calc = EMT()
molecule8.calc = EMT()
e6 = molecule6.get_potential_energy()
e7 = molecule7.get_potential_energy()
e8 = molecule8.get_potential_energy()
print(f'E_6 = {e6} eV, E_7 = {e7} eV, E_8 = {e8} eV')