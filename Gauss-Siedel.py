import numpy as np

def gauss_seidel(A, b, x0, max_iter=100, tol=1e-6):
    """
    Solve the system of linear equations Ax = b using the Gauss-Seidel method.
    
    Parameters:
    A: a matrix of size (n, n)
    b: a vector of size (n,)
    x0: a vector of size (n,) containing the initial guess for the solution
    max_iter: the maximum number of iterations to perform
    tol: the tolerance for the solution (the difference between the current and previous solutions)
    
    Returns:
    x: a vector of size (n,) containing the solution
    """
    x = x0.copy()  # create a copy of the initial guess
    for i in range(max_iter):  # iterate up to the maximum number of iterations
        for j in range(len(x)):  # for each element in the solution
            s = sum(A[j][k] * x[k] for k in range(len(x)) if k != j)  # compute the weighted sum of the neighbors
            x[j] = (b[j] - s) / A[j][j]  # compute the new value for the element
        if np.linalg.norm(x - x0) < tol:  # check for convergence
            return i+1, x
        x0 = x.copy()  # update the initial guess
    return max_iter, x


A = np.array([[4, 1, 1], [1, 3, 1], [1, 1, 2]])
b = np.array([7, 8, 5])
x0 = np.array([1, 1, 1])
its = 0

its, x = gauss_seidel(A, b, x0)

print('Required Iterations:', its)
print('Solution:', x)
