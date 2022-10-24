import numpy as np
import sys



# Conjugate Gradient Solver
def CG_Solver(A,b,x0,tolerance):
    xk = x0
    rk = np.dot(A, xk) - b
    pk = -rk
    rk_norm = np.linalg.norm(rk)
    
    num_iter = 0
    curve_x = [xk]
    while rk_norm > tolerance:
        apk = np.dot(A, pk)
        rkrk = np.dot(rk, rk)
        
        alpha = rkrk / np.dot(pk, apk)
        xk = xk + alpha * pk
        rk = rk + alpha * apk
        beta = np.dot(rk, rk) / rkrk
        pk = -rk + beta * pk
        
        num_iter += 1
        curve_x.append(xk)
        rk_norm = np.linalg.norm(rk)
        
    return np.array(curve_x)

def CG_Solver2(A,b,x0,tolerance):
    rk = b - np.dot(A,x0)
    pk = rk
    rk_norm = np.linalg.norm(rk)
    k = 0
    while rk_norm > tolerance:
        ak = np.dot(rk,rk)/np.dot(pk,np.dot(A,pk))
        xk = xk + ak * pk
        rk2 = rk - ak * np.dot(A,pk)
        rk_norm = np.linalg.norm(rk)
        bk = np.dot(rk2,rk2)/np.dot(rk,rk)
        pk = rk2 + bk * pk
        k = k + 1

    return np.array(xk)
