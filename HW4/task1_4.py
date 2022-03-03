import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('LrTDDFTresults.dat',skiprows=6)
print(data.shape())