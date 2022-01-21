import numpy as np
import matplotlib.pyplot as plt

# HW1 Task 1
def chi(alpha, r):
    return np.exp(-alpha * r**2)

def phi(Cs,alphas, r):
    phi = 0
    for i,alpha in enumerate(alphas):
        phi += Cs[i] * chi(alpha,r)
    return phi
    
r = np.linspace(0,1,10)
print(r)

alphas = [0.297104, 1.236745, 5.749982, 38.216677]