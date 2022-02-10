import ase as ase
from ase.units import fs, kB 

atoms = ase.read('cluster24_9176_Na.xyz')
distances = atoms.get_distances()
for atom in atoms:
    if atom.symbol() == 'Na':
        print('Saltiness detected')
    