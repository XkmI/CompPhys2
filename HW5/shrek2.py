import numpy as np
import matplotlib.pyplot as plt
from ase.io import Trajectory #read
from ase.units import kJ
from ase.eos import EquationOfState

atomNames = ['Au', 'Pt', 'Rh']
latPms = [4.15, 3.95, 3.85]
optInds = [28,27,23]

for i,aN in enumerate(atomNames):
    configs = Trajectory(aN + '.traj') # @0:50 read 50 configurations

    # Extract lattice parameters and energies:
    latPmsVar = np.linspace(0.95*latPms[i],1.05*latPms[i],50)# [np.cbrt(4*atoms.get_volume()) for atoms in configs]
    energies = [atoms.get_potential_energy() for atoms in configs]
    
    plt.plot(latPmsVar,energies,'-', label='Energy for ' + aN + ' lattice')
    plt.plot(latPmsVar[optInds[i]],energies[optInds[i]],'kx',label='Optimal lattice parameter for ' + aN + f', {latPmsVar[optInds[i]]:.2f}')
    #plt.plot(latPms,energies,'kx')
    plt.ylabel('Energy [eV]')
    plt.xlabel('Lattice parameter [Å]')

plt.legend()
plt.savefig('AuPtRh.png')
plt.show()