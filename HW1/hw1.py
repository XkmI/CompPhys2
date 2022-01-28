import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import scipy.integrate as ig
import scipy.sparse as sp

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

def integrand(r2vec, r1vec, Cs, alphas):
    return phi(Cs, alphas, la.norm(r2vec))**2  / la.norm(r1vec-r2vec)

def task1():
    # Define alphas and initial Cs
    alphas = [0.297104, 1.236745, 5.749982, 38.216677]
    Cs = np.ones(4)

    # Initialise and calculate the 4x4 matrices S and h
    SMat = np.zeros([4,4])
    hMat = np.zeros([4,4])
    for p in range(4):
        for q in range(4):
            aSum = alphas[p] + alphas[q]
            aDiv = alphas[p]*alphas[q]/aSum**2
            hMat[p,q] = 4*np.pi * ( 3.0/4.0 * np.sqrt(np.pi / aSum) * aDiv - 1.0/aSum)
            SMat[p,q] = np.power(np.pi/(aSum),1.5)

    # Initialise and calculate the 4x4x4x4 matrix Q (Note that our index ordering is pqrs)
    QMat = np.zeros([4,4,4,4])
    for p in range(4):
        for q in range(4):
            for r in range(4):
                for s in range(4):
                    QMat[p,q,r,s] = 1.0 / ( (alphas[p] + alphas[q])*(alphas[r] + alphas[s])*np.sqrt(alphas[p] + alphas[q]+alphas[r] + alphas[s]) )
    QMat *= 2*np.power(np.pi,2.5)

    # Normalise Cs
    Cnorm2 = np.dot(np.matmul(Cs, SMat), Cs)
    Cs /= np.sqrt(Cnorm2)
    print('------------------')
    print('Initial Cvec is:')
    print(Cs)

    # This will store the ground state energies for the two last runs:
    groundE = np.array([0.0, 1.0])
    i = 0

    while np.abs(groundE[1] - groundE[0]) > 1e-5:
        i += 1
        groundE[0] = groundE[1]

        # Calculate the 4x4 matrix F
        FMat = hMat + np.matmul(np.matmul(Cs,QMat), Cs.transpose())

        # Calculate new Cs
        _, Cs = la.eigh(FMat,SMat, subset_by_index=[0, 0])
        Cs = Cs.transpose()[0]
        print('------------------')
        print('Iteration ' + str(i))
        print('New Cvec is:')
        print(Cs)
        print('Is Cvec normalised? This should be = 1:')
        print(np.dot(np.matmul(Cs, SMat), Cs))

        # Calculate new ground state energy
        groundE[1] = np.dot(np.matmul(Cs,hMat+FMat),Cs)
        print('Ground state energy for two latest iterations are:')
        print(groundE)

'''
xyzMax = 10.0
nPoints = 1001
xyzVals = np.linspace(-xyzMax, xyzMax, nPoints)

vlGrid = np.zeros([nPoints,nPoints,nPoints])
for i in range(nPoints):
    for j in range(nPoints):
        for k in range(nPoints):
            r1vec = np.array([xyzVals[i],xyzVals[j],xyzVals[k]])
            r1 = la.norm(r1vec)

            def currIntegrand(z2,y2,x2):
                r2vec = np.array([x2,y2,z2])
                return integrand(r2vec,r1vec,Cs,alphas)

            vlGrid[i,j,k] = -lapPhi(Cs, alphas, r1)/2.0 + ( - 2.0/r1 + ig.tplquad(currIntegrand,-xyzMax,xyzMax,-xyzMax,xyzMax,-xyzMax,xyzMax) ) * phi(Cs, alphas, r1)
'''
'''
xs = np.linspace(0,10)
plt.plot(xs,phi(Cs, alphas, xs))
plt.show()
'''

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

# HW1 Task 3
def getAMat(vHxc,rs):
    nMin2 = len(rs)-2
    A = sp.diags([-0.5, 1, -0.5], [-1, 0, 1], shape=(nMin2,nMin2)).toarray() / (rs[1]-rs[0])**2

    diagVec = vHxc - 2/rs[1:-1]
    A += np.diag(diagVec)

    return A

def task3():
    rMax = 20
    nPoints = 2001
    rs = np.linspace(0,rMax,nPoints)
    vHxc = 1/rs[1:-1]

    AMat = getAMat(vHxc,rs)

    ei=0
    eps, uVec = la.eigh(AMat,subset_by_index=[ei,ei])
    eps = eps[0]
    uVec = uVec.flatten()
    print(eps)
    print(uVec[0])
    print(uVec[-1])
    plt.plot(rs[1:-1], uVec/(np.sqrt(rs[1]-rs[0])),linewidth=6, label='Calculated u(r)')
    plt.plot(rs[1:-1],2*rs[1:-1]*np.exp(-rs[1:-1]),'-.',linewidth=2, label='Theoretical u(r)')
    plt.xlabel('Radial distance $r$ [$a_0$]')
    plt.ylabel('Radial wave function $u(r)$ [$a_0^{-1/2}$] ')
    plt.xlim([0,10])
    plt.legend()
    plt.show()
    

    

#task1()
#task2()
task3()
