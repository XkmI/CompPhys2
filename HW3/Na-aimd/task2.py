import ase as ase
from ase.units import fs, kB 

def rdf(filename, symbol):
    atoms = ase.read(filename)
    distances = atoms.get_distances()
    for atom in atoms:
        if atom.symbol() == symbol:
            print('Saltiness detected')
    