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

def integrand(z2,y2,x2 , r1vec, Cs, alphas):
    r2vec=np.array([x2,y2,z2])
    return phi(Cs, alphas, la.norm(r2vec))**2  / la.norm(r1vec-r2vec)

alphas = [0.297104, 1.236745, 5.749982, 38.216677]
Cs = np.ones(4)

r1vec = 

-lapPhi(Cs, alphas, r1)/2 + ( - 2/r1 + ig.tplquad( integrand(z2,y2,x2,r1vec,Cs,alphas),-100,100) ) * phi(Cs, alphas, r1)
