import numpy as np

def gen_rand_sym_mat(N):
    """Generates a NxN symmetrical matrix of random real numbers between -10 and 10"""
    M = np.ones((N,N))
    for i in range(N):
        M[i,i] = 20 * np.random.random_sample(size = None) - 10
        for j in range(i):
            M[i,j] = 20 * np.random.random_sample(size = None) - 10
            M[j,i] = M[i,j]
    return M

def gen_sym_mat(N):
    """Generates a NxN matrix of random real numbers between -10 and 10"""
    M = np.ones((N,N))
    for i in range(N):
        for j in range(N):
            M[i,j] = 20 * np.random.random_sample(size = None) - 10
    return M

def gen_rand_tridiag(N):
    """Generates a random NxN tridiagonal matrix of real numbers between -10 and 10"""
    M = np.zeros((N,N))
    for i in range(N):
        if i == 0:
            for j in range(2):
                M[i,j] = 20 * np.random.random_sample(size = None) - 10
        elif i == N-1:
            for j in range(i-1,i+1):
                M[i,j] = 20 * np.random.random_sample(size = None) - 10
        else:
            for j in range(i-1,i+2):
                M[i,j] = 20 * np.random.random_sample(size = None) - 10
    print(M)
    return M
