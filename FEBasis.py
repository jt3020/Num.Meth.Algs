import numpy as np
import matplotlib.pyplot as plt


#Finite Element basis functions and derivatives
def Eval_FE(p,x):
    N_Xi = np.zeros(p+1)
    XiVec = np.zeros(p+1)
    for i in range (p+1):
        XiVec[i] = (1.0/p)*float(i)
    
    for j in range(p+1):
        N_Xi[j] = 1.0
        for i in range(p+1):
            if i != j:
                N_Xi[j] = N_Xi[j]*(x-XiVec[i])/(XiVec[j]-XiVec[i])

    return N_Xi


def plot_splines(degree, nx):
    
    # Grid points for evaluation
    xs = np.linspace(0.0,1.0,nx)

    # Number of the BSplines
    N = degree + 1

    ys = np.zeros((N,nx))
    for ix,x in enumerate(xs):
        N_Xi    = Eval_FE(degree,x)
        ys[:, ix] = N_Xi[:]
        
    for i in range(0,N):
        plt.plot(xs,ys[i,:], label='$N_{}$'.format(i+1))
    plt.legend(loc=9, ncol=4)
    plt.xlabel("x")
    plt.ylabel("N")
    plt.show()

degree = 2
nx = 100
plot_splines(degree, nx)




