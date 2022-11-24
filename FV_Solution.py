import numpy as np
import matplotlib.pyplot as plt

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

# Problem Size
n = 6
Delta = 1/(n-1)

# Material Properties
D = 1.0/3.0
Sig_a = 1.0
S = 1.0

# Matrices and Vectors associated with problem
A = np.zeros((n,n))
b = np.zeros(n)
x = np.zeros(n)

# Fill in problem
# First Node
A[0,1] = -D/Delta
A[0,0] = Sig_a*Delta - A[0,1]
b[0] = S*Delta
# Internal Nodes
for i in range(1,n-1):
    A[i,i-1] = -D/Delta
    A[i,i+1] = -D/Delta
    A[i,i] = Sig_a*Delta - (A[i,i-1]+A[i,i+1])
    b[i] = S*Delta
# Last Node
A[n-1,n-2] = -D/Delta
A[n-1,n-1] = Sig_a*Delta - A[n-1,n-2]
b[n-1] = S*Delta

# Apply Boundary Conditions
A[0,0] = 1E5
A[n-1,n-1] = 1E5

# Solve the problem
x = T_alg(A,b,x)

# Generate flux plot
position = np.zeros(n)
for i in range(n):
    position[i] = i*Delta
plt.plot(position,x)
plt.xlabel("Position")
plt.ylabel("Flux")
plt.show()
