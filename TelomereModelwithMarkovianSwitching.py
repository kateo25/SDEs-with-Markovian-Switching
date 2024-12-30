import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from numpy.random import random_sample
from scipy.stats import rv_discrete

a1 = 0.41*10**(-6)
c1 = 7.5

a2 = 0.22*10**(-6)
c2 = 4.5

def drift (s, x):
    if s == 0:
        f = -c1 - (0.5*a1*x**2)
    elif s == 1:
        f = -c2 - (0.5*a1*x**2)
    elif s == 2:
        f = -c1 - (0.5*a2*x**2)
    else:
        f = -c2 - (0.5*a2*x**2)
    return f

def diffusion (s, x):
    if s == 0:
        g = np.sqrt((a1*x**3)/3)
    elif s == 1:
        g = np.sqrt((a1*x**3)/3)
    elif s == 2:
        g = np.sqrt((a2*x**3)/3)
    else:
        g = np.sqrt((a2*x**3)/3)
    return g

def diffusionPrime (s, x):
    if s == 0:
        gprime = (np.sqrt(3*a1*x**3))/(2*x)
    if s == 1:
        gprime = (np.sqrt(3*a1*x**3))/(2*x)
    if s == 2:
        gprime = (np.sqrt(3*a2*x**3))/(2*x)
    else:
        gprime = (np.sqrt(3*a2*x**3))/(2*x)
    return gprime


s = (0,1,2,3)
A = np.array([[-0.3, 0.1, 0.1, 0.1], [0.1, -0.3, 0.1, 0.1], [0.1, 0.1, -0.3, 0.1], [0.1, 0.1, 0.1, -0.3]])
print(A)
seed = 100
rng=np.random.default_rng(seed)
B=100
Xend = np.array([])


for z in range(B):
    
    i = 0
    j = 0
    k = 0
    l = 0
    q = 0
    tee= 0
    t_comb = np.array([0])
    s_comb = np.array([0])
    ess = int(s_comb[0])
    Xtemp = 1000
    Xem = np.array([])
    tvec = np.array([0])

    while t_comb[q] < 31:

        lambdai =  -A[i,i]
        b = 1/lambdai
        pi = np.random.default_rng().exponential(scale=b)
        tee = tee + pi
        t_comb = np.append(t_comb, [tee])
        q = q+1
        p = np.empty(0)

        for y in range(len(s)):
            if y != i:
                x = A[i, y]/ -A[i,i]
                p = np.append(p, x)
            else:
                pass
                
        pList = list(p)  
        svec = np.empty(0)

        for y in range(len(s)):
            if y < i:
                svec = np.append(svec, y)
            elif y > i:
                svec = np.append(svec, y)
            else:
                pass

        values = svec
        distrib = rv_discrete(values=(range(len(values)), pList))  
        distrib.rvs(size=1) 
        a = values[_]
        state = a[0]
        s_comb = np.append(s_comb, [state])
        i = int(state)

    tau = np.delete(t_comb, [0])
    tau = np.sort(np.append(tau, [30]))

    while tvec[j] < 30: 

        hmin = 0.002
        hmax = 0.03
        kappa = 10
        v = hmax/(Xtemp**(1/kappa))
        taucheck = tau[l] - tvec[j]

        if v < hmax and hmin < v and v < taucheck:
            hnplusone = v
            tvec = np.append(tvec, [hnplusone + tvec[j]])
        elif hmax < v and hmin < hmax and hmax < taucheck:
            hnplusone = hmax
            tvec = np.append(tvec, [hnplusone + tvec[j]])
        elif v < hmax and hmin > v and hmin < taucheck:
            hnplusone = hmin
            tvec = np.append(tvec, [hnplusone + tvec[j]])
        elif hmax < v and hmin > hmax and hmin < taucheck:
            hnplusone = hmin
            tvec = np.append(tvec, [hnplusone + tvec[j]])
        elif taucheck < v and taucheck < hmax and taucheck < hmin:
            hnplusone = taucheck
            tvec = np.append(tvec, [hnplusone + tvec[j]])
            l = l+1
        elif taucheck < v and taucheck < hmax and taucheck > hmin:
            hnplusone = taucheck
            tvec = np.append(tvec, [hnplusone + tvec[j]])
            l = l+1
        else:  
            print("no, maybe hmin greater than hmax, otherwise idk")

        DW = np.sqrt(hnplusone)*rng.normal(0,1,1)
        if float(tvec[j]) == float(tau[k]):
            if taucheck < hmin:
                print('backstop')
                ess = int(s_comb[k+1])
                k = k+1
            else:
                Xtemp = Xtemp + hnplusone*drift(ess, Xtemp) + diffusion(ess, Xtemp)*DW + (1/2)*diffusionPrime(ess, Xtemp)*diffusion(ess, Xtemp)*(DW**2 - hnplusone)
                ess = int(s_comb[k+1])
                k = k+1
        else:
            Xtemp = Xtemp + hnplusone*drift(ess, Xtemp) + diffusion(ess, Xtemp)*DW + (1/2)*diffusionPrime(ess, Xtemp)*diffusion(ess, Xtemp)*(DW**2 - hnplusone)
            
        Xem = np.append(Xem, [Xtemp]) 
        j = j+1     
    Xend = np.append(Xend, Xem[-1])


Xend
tvec
tau
Xem[-1]
len(Xem)
len(Xend)
len(tvec)
len(tau)



hmin = 0.002
hmax = 0.03
kappa = 10
v = hmax/(Xtemp**(1/kappa))
v

np.mean(Xend)


fig, ax1 = plt.subplot(figsize=(15,10))
#plt.grid(axis='both', color='0.95')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.hist(Xend, bins = 25, ec = 'black', density = True)
plt.xlabel("Telomere Length", fontsize = 15)
plt.ylabel("Density", fontsize=15)
plt.show()



