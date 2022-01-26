import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import scipy.integrate as ig

# HW1 Task 1
def chi(alpha, r):
    return np.exp(-alpha * r**2)

def phi(Cs, alphas, r):
    phi = 0
    for i,alpha in enumerate(alphas):
        phi += Cs[i] * chi(alpha,r)
    return phi

def lapPhi(Cs, alphas, r):
    lapPhi = 0
    for i,alpha in enumerate(alphas):
        lapPhi += 2 * alpha * (2 * alpha * r**2 - 3) * Cs[i] * chi(alpha,r)
    return lapPhi 

def H(Cs, alphas, rvec):
   -lapPhi(Cs, alphas, r)/2 - 2/r + ig.quad(np.abs(phi(Cs, alphas, np.abs(rvec)))**2 / np.abs(rvec-Rvec))


r = np.linspace(0,1,10)
print(r)

alphas = [0.297104, 1.236745, 5.749982, 38.216677]
Cs = np.ones(4)
