import ase.io as io
import numpy as np
from ase.units import fs, kB 

def partial_rdf(filename):
    atoms = io.read(filename)
    n = 100
    distances = np.zeros(24) #24 oxygen atoms (and water molecules)
    for i in range(49,72): #Oxygen, this is pretty much the center of each water molecule
        distances[i-48] = atoms.get_distances(i,48)
    distances.sort()

    r = np.linspace(0,max(distances),n)  

        

partial_rdf('cluster24_9176_Na.xyz')
