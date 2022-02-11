from ase.io import *
from ase.io.trajectory import TrajectoryReader
import os

# This code is partly stolen from https://sites.lsa.umich.edu/zimmerman-lab/tutorial/...
# ...surface-growing-string-method/converting-trajectory-file-to-xyz-format/

# Read in .traj file
traj = TrajectoryReader("NaCluster24.traj")

# Where to start and how many snapshots
startInd = 9000
nSnapshots = 1

path = os.getcwd()
#path = path + '/scratch'
#if not os.path.exists(path): os.makedirs(path)

# Output file name
outFileName = 'NaCluster24.xyz'
# Write each selected structure from the .traj file in .xyz format
for i in range(startInd, startInd+nSnapshots):
    atoms = traj[i]
    string = 'structure%03d' % (i,) +'.xyz'
    outStruct = os.path.join(path, string)
    write(outStruct, atoms)
    
    # Combine all structures in one .traj file
    inFile = open(os.path.join(path, 'structure%03d' % (i,)  +'.xyz'), 'r')
    fileStr = inFile.read()
    outFile = open(outFileName, 'a')
    outFile.write(fileStr)