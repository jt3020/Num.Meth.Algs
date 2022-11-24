import numpy as np
import sys
import os
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt


def find_span( knots, degree, x ):
    # Knot index at left/right boundary
    low  = degree
    high = 0
    high = len(knots)-1-degree

    # Check if point is exactly on left/right boundary, or outside domain
    if x <= knots[low ]: returnVal = low
    elif x >= knots[high]: returnVal = high-1
    else:
        # Perform binary search
        span = (low+high)//2
        while x < knots[span] or x >= knots[span+1]:
            if x < knots[span]:
                high = span
            else:
                low  = span
            span = (low+high)//2
        returnVal = span

    return returnVal


def all_bsplines( knots, degree, x, span ):
    left   = np.empty( degree  , dtype=float )
    right  = np.empty( degree  , dtype=float )
    values = np.empty( degree+1, dtype=float )

    values[0] = 1.0
    for j in range(0,degree):
        left [j] = x - knots[span-j]
        right[j] = knots[span+1+j] - x
        saved    = 0.0
        for r in range(0,j+1):
            temp      = values[r] / (right[r] + left[j-r])
            values[r] = saved + right[r] * temp
            saved     = left[j-r] * temp
        values[j+1] = saved

    return values


def plot_splines(knots, degree, nx=100):
    xmin = knots[degree]
    xmax = knots[-degree-1]
    
    # grid points for evaluation
    xs = np.linspace(xmin,xmax,nx)

    # this is the number of the BSplines in the Schoenberg space
    N = len(knots) - degree - 1

    ys = np.zeros((N,nx), dtype=np.double)
    for ix,x in enumerate(xs):
        span = find_span( knots, degree, x )    
        b    = all_bsplines( knots, degree, x, span )  
        ys[span-degree:span+1, ix] = b[:]
        
    for i in range(0,N):
        plt.plot(xs,ys[i,:], label='$N_{}$'.format(i+1))
    plt.legend(loc=9, ncol=4)
    plt.show()


def point_on_bspline_curve(knots, P, x):
    degree = len(knots) - len(P) - 1
    d = P.shape[-1]

    span = find_span( knots, degree, x )
    b    = all_bsplines( knots, degree, x, span )

    c = np.zeros(d)
    for k in range(0, degree+1):
        c[:] += b[k]*P[span-degree+k,:]
    return c   


def Curve_Plot(knots,degree,C,nx):
    n      = len(knots) - degree - 1  

    xs = np.linspace(0., 1., nx)
    
    Q = np.zeros((nx, 2))
    for i,x in enumerate(xs):
        Q[i,:] = point_on_bspline_curve(knots, C, x)

    plt.plot(Q[:,0], Q[:,1], '-b')
    plt.plot(C[:,0], C[:,1], '--ok', linewidth=0.7)
    
    for i in range(0, n):
        x,y = C[i,:]
        plt.text(x+0.05,y+0.05,'$\mathbf{C}_{' + str(i) + '}$')

    plt.axis([-2.2, 1.2, -1.2, 1.2])

    plt.show()


knots  = [0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1.]
degree = 2
n      = len(knots) - degree - 1

C = np.zeros((n, 2))
C[:, 0] = [0.,  1., 1., -2., -1., -2.]
C[:, 1] = [.5, -1., 1.,  1.,  0., -1.]

Curve_Plot(knots,degree,C,1000)

# T = np.array([0, 0, 0, 1, 2, 3, 4, 5, 5, 5])

# plot_splines(T, degree=2, nx=100)
# plt.show()




