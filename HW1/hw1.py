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


# HW1 Task 2
def lap(U,i,h):
    return (U[i+1]-2*U[i] + U[i-1])/h**2

def realphi(r):
    return (1/r) - (1 + 1/r) * np.exp(-2*r)

def task2():
    a=0 
    b=10 
    n=1000 
    h=(b-a)/n
    rmax=b #idfk

    r = np.linspace(a,b,n)

    U = np.zeros(n)
    ddU = np.zeros(n)
    phisquare = np.zeros(n)

    for i in range(n):
        if i == 0 or n:
            U[i] = 0
        U[i]=r[i] * V[i] - r[i]/rmax #detta är sus
    
    for i in range(n):
        if i == 0 or n:
            U[i] = 0
        ddU[i] = lap(U,i,h)

    phisquare = - 1/(4*np.pi) * ddU

    plt.plot(r,np.sqrt(phisquare),'r')
    plt.plot(r,realphi(r),'b')
    plt.show()

task1()
task2()
task3()
task4()
task5()
task6()