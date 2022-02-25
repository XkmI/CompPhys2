from ase.io import read
from ase.units import kJ
from ase.utils.eos import EquationOfState

atomNames = ['Au', 'Pt', 'Rh']
latPms = [4.05, 3.90, 3.80]

for aN in atomNames:
    configs = read(aN + '.traj@0:7') # read 7 configurations

    # Extract volumes and energies:
    volumes = [atoms.get_volume() for atoms in configs]
    energies = [atoms.get_potential_energy() for atoms in configs]
    eos = EquationOfState(volumes , energies)
    v0, e0, B = eos.fit()

    print(e0)
    print(f'{v0} \AA^3, {B / kJ * 1.0e24} GPa')
    eos.plot(aN + '-eos.png')