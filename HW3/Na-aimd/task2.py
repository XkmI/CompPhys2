import ase.io as io
import numpy as np
from ase.units import fs, kB 
import matplotlib.pyplot as plt
from ase.io.trajectory import TrajectoryReader
import os

def partial_rdf(filename=None, atoms=None, dr=None, nBins=None):
    if atoms is None:
        atoms = io.read(filename)
    nOx = 24 #24 oxygen atoms (and water molecules)
    distances = np.zeros(nOx) 
    for i in range(nOx): #Finds distances from the sodium atom (#-1 i.e. the last one) to all oxygen atoms (#0-23), since they are more or less at the centre of each water molecule
        distances[i] = atoms.get_distances(i,-1,mic=True)
    distances.sort()

    # Set dr and #bins for the ~histogram of oxygen distances
    if dr is None:
        dr = 1.0
    if nBins is None:
        nBins = 20 #int(np.ceil(max(distances))) + 1
    # Set limits for histogram bins and calculate centre points of bins
    rLimBins = np.linspace(0,nBins*dr,nBins+1)
    rBins = rLimBins[1:]-0.5*dr

    # Initialise and fill bins appropriately
    nOxHist = np.zeros(nBins)
    for bInd in np.floor(distances/dr):
        nOxHist[int(bInd)] += 1

    prdf = nOxHist * max(distances)**3 / (dr * 3 * nOx * rBins**2)
    #print(max(distances))
    #plt.plot(rBins,prdf)
    #plt.show()

    return (rBins,prdf,nOxHist)

# Read in .traj file
traj = TrajectoryReader("babys_first_NaCluster24.traj") #babys_first_
# Where to start and how many snapshots
startInd = 1000
nSnapshots = 2340
path = os.getcwd()

dr = 0.1
nBins = 200
rLimBins = np.linspace(0,nBins*dr,nBins+1)
rBins = rLimBins[1:]-0.5*dr
rdf = np.zeros(nBins)
nOxHist = np.zeros(nBins)

for i in range(startInd, startInd+nSnapshots):
    atoms = traj[i]
    #print(atoms[0])
    _, prdf, nOxHist_part = partial_rdf(atoms=atoms,dr=dr,nBins=nBins)
    rdf += prdf
    nOxHist += nOxHist_part

rdf /= nSnapshots
nOxHist /= nSnapshots

print(sum(nOxHist[:int(3.2/dr)]))
print(nOxHist[:int(3.2/dr)])

plt.plot(rBins,rdf,label='Radial distribution function of O around Na')

plt.plot([3.15, 3.15],[-5, 20],'k-.', label='First minimum ($r = 3.15 a_0$)')
plt.xlabel('Radial distance $r$ [$a_0$]')
plt.ylabel('RDF$_\\operatorname{NaO}$ at distance r')
plt.xlim([-0.1, 8])
plt.ylim([-0.1, 8])
plt.legend()
plt.show()
plt.savefig('lamao.pdf')


