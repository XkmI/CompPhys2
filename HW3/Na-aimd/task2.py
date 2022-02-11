import ase as ase
from ase.units import fs, kB 

def partial_rdf(filename):
    atoms = ase.read(filename)
    #find the Na atom somehow
    distances = atoms.get_distances()
    for atom in atoms:
        if atom.symbol() == 'H':
            print('Saltiness detected')
        if atom.symbol() == 'O':