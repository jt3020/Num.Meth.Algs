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

#Gaussian Quadrature
def G_Quad(n,weight,absicca):
    if n == 2:
        weight[0] = 0.5
        weight[1] = 0.5
        absicca[0] = 0.211
        absicca[1] = 0.789
    elif n == 3:
        weight[0] = 0.278
        weight[1] = 0.444
        weight[2] = 0.278
        absicca[0] = 0.113
        absicca[1] = 0.5
        absicca[2] = 0.887
    else: 
        print("Error: G_Quad only set up for n < 4")
        exit()

#Finite Element basis functions and derivatives
def Eval_FE(p,Pos,N_Xi,dN_Xi):

    XiVec = np.zeros(p+1)
    for i in range (p+1):
        XiVec[i] = (1.0/p)*i
    
    for j in range(p+1):
        N_Xi[j] = 1.0
        for i in range(p+1):
            if i != j:
                N_Xi[j] = N_Xi[j]*(Pos-XiVec[i])/(XiVec[j]-XiVec[i])

    for j in range(p+1):
        dN_Xi[j] = 0.0
        for i in range(p+1):
            if i != j:
                TempValue = 1.0/(XiVec[j]-XiVec[i])
                for k in range(p+1):
                    if (k != j) and (k != i):
                        TempValue = TempValue*(Pos-XiVec[k])/(XiVec[j]-XiVec[k])
                dN_Xi[j] = dN_Xi[j] + TempValue

#Jacobian Calculation
def Jacobian(dN_Xi,Nodes):
    Det_J = 0.0
    for i in range(2):
        Det_J = Det_J + Nodes[i]*dN_Xi[i]
    Det_J = abs(Det_J)
    return Det_J
    

# Problem Geometry
n = 4
Nodes = np.zeros((n,2))

# Get positions of each elements nodes 
# Element 1
Nodes[0,0] = 0.0
Nodes[0,1] = 0.25
# Element 2
Nodes[1,0] = 0.25
Nodes[1,1] = 0.5
# Element 3
Nodes[2,0] = 0.5
Nodes[2,1] = 0.75
# Element 4
Nodes[3,0] = 0.75
Nodes[3,1] = 1.0


# Material Properties
D = 1.0/3.0
Sig_a = 1.0
S = 1.0

# Matrices and Vectors associated with problem
M = np.zeros((n,n))
K = np.zeros((n,n))
F = np.zeros(n)
x = np.zeros(n)

# Gaussian Quadrature
N_QuadPts = 2
weight = np.zeros(N_QuadPts)
absicca = np.zeros(N_QuadPts)
G_Quad(N_QuadPts,weight,absicca)

# Basis functions and derivatives
N_Xi = np.zeros(2)
dN_Xi = np.zeros(2)

# Loop over elements
for i in range(n):
    # Loop over quadrature points
    for j in range(N_QuadPts):
        Eval_FE(1,absicca[j],N_Xi,dN_Xi)
        Det_J = Jacobian(dN_Xi,Nodes[i,:])

        for ii in range(i,i+1):
            F[ii] = F[ii] + Det_J*N_Xi[ii]*weight[j]
            for jj in range(i,i+1):
                M[ii,jj] = M[ii,jj] + Det_J*N_Xi[ii]*N_Xi[jj]*weight[j]
                K[ii,jj] = K[ii,jj] + Det_J*dN_Xi[ii]*dN_Xi[jj]*weight[j]

# Solve the problem
x = T_alg(M+K,F,x)

# Generate flux plot
position = np.zeros(n)
for i in range(n):
    position[i] = i*0.25
plt.plot(position,x)
plt.xlabel("Position")
plt.ylabel("Flux")
plt.show()
