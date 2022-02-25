import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState

atomNames = ['Au', 'Pt', 'Rh']
latPms = [4.05, 3.90, 3.80]

for aN in atomNames:
    configs = read(aN + '.traj') # @0:15 read 15 configurations

    # Extract lattice parameters and energies:
    latPms = [np.cbrt(4*atoms.get_volume()) for atoms in configs]
    energies = [atoms.get_potential_energy() for atoms in configs]
    
    plt.plot(latPms,energies,'-', label='Energy for ' + aN + ' lattice')
    plt.plot(latPms,energies,'kx')
    plt.ylabel('Energy [eV]')
    plt.xlabel('Lattice parameter [Å]')

plt.legend()
plt.savefig('AuPtRh.png')
plt.show()