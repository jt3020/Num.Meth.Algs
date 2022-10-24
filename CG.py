import numpy as np
import sys

# Conjugate Gradient Solver
def CG_Solver(A,b,x,tolerance):
    r = b - np.dot(A,x)
    p = r 
    rsold = np.dot(r,r)
    k = 0

    while np.linalg.norm(rsold) > tolerance:
        Ap = np.dot(A,p)
        alpha = rsold/np.dot(p,Ap)
        x = x + np.dot(alpha,p)
        r = r - np.dot(alpha,Ap)
        rsnew = np.dot(r,r)
        p = r + np.dot((rsnew/rsold),p)
        rsold = rsnew
        k = k + 1
    return np.array(x)

# Pre-Conditioned Conjugate Gradient Method
def PCG_Solver(A,M,b,x,tolerance):
    r = b - np.dot(A,x)
    z = np.dot(M,r)
    p = z 
    r2 = np.dot(r,r)
    k = 0

    while np.linalg.norm(r2) > tolerance:
        rz = np.dot(r,z)
        Ap = np.dot(A,p)
        alpha = rz/np.dot(p,Ap)
        x = x + np.dot(alpha,p)
        r2 = r - np.dot(alpha,Ap)
        
        z2 = np.dot(M,r2)
        beta = np.dot(r2,z2)/np.dot(r,z)
        p = z2 + np.dot(beta,p)
        k = k + 1

    return np.array(x)

# Bi-Conjugate Gradient Solver
def BCG_Solver(A,M,b,x,tolerance):
    r = b - np.dot(A,x)
    rho1 = np.dot(r,r)
    rh = r
    k = 0

    while np.linalg.norm(rho1) > tolerance:
        if k != 0:
            beta = (rho1/rho2)*(alpha/omega)
            p = r + beta*(p - (omega*v))
        else:
            p = r

        ph = np.dot(M,p)
        v = np.dot(A,ph)

        alpha = rho1/np.dot(rh,v)
        s = r - (alpha*v)

        sh = np.dot(M,s)
        t = np.dot(A,sh) 

        if np.linalg.norm(t) < tolerance:
            omega = 0.0
        else:
            omega = np.dot(t,s)/np.dot(t,t)
        
        x = x + (alpha*ph) + (omega*sh)
        r = s - (omega*t)

        rho1 = np.dot(r,r)
        rho2 = rho1  

        k = k + 1     

    return(np.array(x))


N = 5
x = np.ones(N)
b = np.ones(N)
A = np.zeros((N,N))
M = np.zeros((N,N))

for i in range(5):
    A[i,i] = 2.0
    M[i,i] = 1.0/A[i,i]

x = PCG_Solver(A,M,b,x,1E-5)

print('x:',x)
print('A:',A)
