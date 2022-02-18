from ase import Atoms
from gpaw import GPAW, FermiDirac
import ase.io as io

atoms = io.read('8.xyz')

calc = GPAW(nbands=10, h=0.25, txt='out.txt', occupations=FermiDirac(0.05), setups={'Na': '1'}, mode='lcao', basis='dzp')
atoms.calc = calc

# Start a calculation:
energy = atoms.get_potential_energy()

# Save wave functions:
calc.write('na8.gpw', mode='all')
