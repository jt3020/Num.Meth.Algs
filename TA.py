import numpy as np

#Thomas Algorithm
def T_alg(A,b,x):
    """Solves a tridiagonal system using the Thomas Algorithm"""
    N = A.shape[0]
    for i in range(1,N):
        w = A[i,i-1] / A[i-1,i-1]
        A[i,:] = A[i,:] - A[i-1,:] * w
        b[i] = b[i] - b[i-1] * w
    for i in range(N-1,-1,-1):
        if i == N-1:
            x[i] = b[i] / A[i,i]
        else:
            x[i] = (b[i] - A[i,i+1] * x[i+1]) / A[i,i]
    return x
