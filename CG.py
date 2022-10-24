import numpy as np
import sys

# Conjugate Gradient Solver
def CG_Solver(A,b,x,tolerance):
    r = b - np.dot(A,x)
    p = r 
    rsold = np.dot(r,r)

    while np.linalg.norm(rsold) > tolerance:
        Ap = np.dot(A,p)
        alpha = rsold/np.dot(p,Ap)
        x = x + np.dot(alpha,p)
        r = r - np.dot(alpha,Ap)
        rsnew = np.dot(r,r)
        p = r + np.dot((rsnew/rsold),p)
        rsold = rsnew
    return np.array(x)


N = 5
x = np.ones(N)
b = np.ones(N)
A = np.zeros((N,N))

for i in range(5):
    A[i,i] = 2.0

x = CG_Solver(A,b,x,1E-5)

print('x:',x)
print('A:',A)
