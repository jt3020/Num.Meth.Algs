import numpy as np

def eigenvalue_solver(A):
    """
    Find the eigenvalues and eigenvectors of the matrix A.

    Parameters:
    A: a square matrix of size (n, n)

    Returns:
    eigenvalues: a list of size (n,) containing the eigenvalues of A
    eigenvectors: a list of size (n,) containing the eigenvectors of A
    """
    n = len(A)
    eigenvalues = []
    eigenvectors = []
    for i in range(n):  # for each row in the matrix
        # solve the equation (A - lambda*I)x = 0 for x
        matrix = A - np.eye(n) * A[i][i]
        vector = np.zeros(n)
        x = np.linalg.solve(matrix, vector)
        eigenvalues.append(A[i][i])
        eigenvectors.append(x)
    return eigenvalues, eigenvectors

A = [[4, 1, 1], [1, 3, 1], [1, 1, 2]]
eigenvalues, eigenvectors = eigenvalue_solver(A)

print('Matrix A:',A)
print('Eigenvalues:',eigenvalues)
print('Eigenvectors:',eigenvectors)
