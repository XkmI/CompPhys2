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
    
    rs = np.linspace(0,8,1000)
    
    plt.plot(rs, phi(Cs,alphas,rs), label='Optimal wave function $\phi(r)$')
    plt.legend()
    plt.xlabel('Radial distance $r$ [$a_0$]')
    plt.ylabel('3D probability density [$a_0^{-3/2}$]')
    '''
    plt.plot(rs, 2*np.sqrt(np.pi) * rs * phi(Cs,alphas,rs),label='Optimal radial wave function $2\sqrt{\pi} r \phi(r)$')
    plt.legend()
    plt.xlabel('Radial distance $r$ [$a_0$]')
    plt.ylabel('1D probability density [$a_0^{-1/2}$]')
    '''
    plt.show()


# HW1 Task 2
def getVH(r):
    return (1/r) - (1 + 1/r) * np.exp(-2*r)

def task2():
    rMax = 20
    nPoints = 2001
    rs = np.linspace(0,rMax,nPoints)
    nMin2 = nPoints-2
    A = sp.diags([1, -2, 1], [-1, 0, 1], shape=(nMin2,nMin2)).toarray() / (rs[1] - rs[0])**2
    #print(A)
    
    f = - 4 * rs[1:-1] * np.exp(-2*rs[1:-1])
    
    U0vec = la.solve(A,f)
    
    VHvec = U0vec/rs[1:-1] + 1/rMax

    plt.plot(rs[1:-1],VHvec,linewidth=6,label="Calculated Hartree potential")
    plt.plot(rs,getVH(rs),'-.', label="Theoretical Hartree potential")
    plt.legend()
    plt.xlabel('Radial distance $r$ [$a_0$]')
    plt.ylabel('Potential energy $V_H(r)$ [eV]')
    plt.show()
    #plt.savefig('task2.pdf')

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
    

# HW1 task 4
def task4():
    rMax = 20
    nPoints = 4001
    nMin2 = nPoints-2
    rs = np.linspace(0,rMax,nPoints)
    epsList = np.array([0.0, 1.0])
    i = 0

    # initial condition
    u = 2 * rs[1:-1] * np.exp(-rs[1:-1]) * np.sqrt(2)

    while np.abs(epsList[1] - epsList[0]) > 1e-5:
        i+=1
        epsList[0] = epsList[1]
        print('------------------')
        print('Iteration ' + str(i))

        f = -u**2/rs[1:-1]
        Amat = sp.diags([1, -2, 1], [-1, 0, 1], shape=(nMin2,nMin2)).toarray() / (rs[1] - rs[0])**2
        U0vec = la.solve(Amat,f)
        vsH = U0vec/rs[1:-1] + 1/rMax  #*2?
        BMat = getAMat(vsH,rs)

        ei=0
        eps, u = la.eigh(BMat,subset_by_index=[ei,ei])
        epsList[1] = eps[0]
        u = u.flatten()/(np.sqrt(rs[1]-rs[0]))

        #print(np.sum(u**2))
        print('Eigenvalues for two latest iterations are:')
        print(epsList)
    
    print('------------------')
    print('Calculated ground state energy:')
    eG = 2*eps[-1] - np.sum(u**2 * vsH)*(rs[1]-rs[0])
    print(eG)
    print('------------------')
    print('u next to endpoints:')
    print(u[0])
    print(u[-1])
    plt.plot(rs[1:-1], u, label='Calculated u(r)') # ,linewidth=6
    #plt.plot(rs[1:-1],2*rs[1:-1]*np.exp(-rs[1:-1]),'-.',linewidth=2, label='Theoretical u(r)')
    plt.xlabel('Radial distance $r$ [$a_0$]')
    plt.ylabel('Radial wave function $u(r)$ [$a_0^{-1/2}$] ')
    plt.xlim([0,8])
    plt.legend()
    plt.show()
    

#task1()
#task2()
#task3()
task4()