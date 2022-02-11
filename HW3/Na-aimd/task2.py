import ase.io as io
import numpy as np
from ase.units import fs, kB 
import matplotlib.pyplot as plt

def partial_rdf(filename):
    atoms = io.read(filename)
    nOx = 24 #24 oxygen atoms (and water molecules)
    distances = np.zeros(nOx) 
    for i in range(nOx): #Finds distances from the sodium atom (#48) to all oxygen atoms, since they are more or less at the centre of each water molecule
        distances[i] = atoms.get_distances(49+i,48)
    distances.sort()

    # Choose dr and #bins for ~histogram of oxygen distances
    dr = 1.0
    nBins = 15 #int(np.ceil(max(distances))) + 1
    # Set limits for histogram bins and calculate centre points of bins
    rLimBins = np.linspace(0,nBins*dr,nBins+1)
    rBins = rLimBins[1:]-0.5*dr

    # Initialise and fill bins appropriately
    nOxHist = np.zeros(nBins)
    for bInd in np.floor(distances/dr):
        nOxHist[int(bInd)] += 1

    prdf = nOxHist * max(distances)**3 / (dr * 3 * nOx * rBins**2)
    
    #plt.plot(rBins,prdf)
    #plt.show()

    return prdf
    


