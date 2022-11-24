import numpy as np
import matplotlib.pyplot as plt

def find_span(knots, degree, x):
    # Knot index at left/right boundary
    low  = degree
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

def all_bsplines(knots, degree, x, span):
    left   = np.zeros(degree)
    right  = np.zeros(degree)
    values = np.zeros(degree+1)

    values[0] = 1.0
    for j in range(degree):
        left [j] = x - knots[span-j]
        right[j] = knots[span+1+j] - x
        saved    = 0.0
        for r in range(j+1):
            temp      = values[r] / (right[r] + left[j-r])
            values[r] = saved + right[r] * temp
            saved     = left[j-r] * temp
        values[j+1] = saved

    return values

def plot_splines(knots, nx=100):
    # Get Degree of knot vector
    for i in range(knots.shape[0]):
        if knots[i] > 0:
            degree = i-1
            break

    xmin = knots[degree]
    xmax = knots[-degree-1]
    
    # Grid points for evaluation
    xs = np.linspace(xmin,xmax,nx)

    # Number of the BSplines
    N = len(knots) - degree - 1

    ys = np.zeros((N,nx), dtype=np.double)
    for ix,x in enumerate(xs):
        span = find_span(knots, degree, x)    
        b    = all_bsplines(knots, degree, x, span)  
        ys[span-degree:span+1, ix] = b[:]
        
    for i in range(0,N):
        plt.plot(xs,ys[i,:], label='$N_{}$'.format(i+1))
    plt.legend(loc=9, ncol=4)
    plt.xlabel("x")
    plt.ylabel("N")
    plt.show()

knots = np.array([0, 0, 0, 0.2, 0.5, 0.5, 0.8, 1, 1, 1])

plot_splines(knots, nx=100)




