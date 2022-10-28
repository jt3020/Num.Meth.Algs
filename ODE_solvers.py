import numpy as np

def Backward_Euler(func,y0,t0,tf,tolerance,delta):
    """Calculate the solution of func(y,t) at tf with initial condition y0 at t0"""
    delta_t = tf - t0
    yf = y0
    residual = yf - y0 - delta_t * func(yf,tf)
    while np.absolute(residual) > tolerance:
        
        m = ((yf + delta - y0 - delta_t * func(yf + delta,tf)) - (yf - y0 - delta_t * func(yf,tf))) / delta
        yf = yf - (yf - y0 - delta_t * func(yf,tf)) / m
        residual = yf - y0 - delta_t * func(yf,tf)
        print(residual)
    return yf


def RK45(func,y0,t0,tf,tolerance):
    """Calculate the solution of func(y,t) at tf with initial condition y0 at t0"""
    delta_t = tf - t0
    a = np.array([[0,0,0,0,0,0],
                  [1/4,0,0,0,0,0],
                  [3/32,9/32,0,0,0,0],
                  [1932/2197,-7200/2197,7296/2197,0,0,0],
                  [439/216,-8,3680/513,-845/4104,0,0],
                  [-8/27,2,-3544/2565,1859/4104,-11/40,0]])
    b5 = np.array([16/135,0,6656/12825,28561/56430,-9/50,2/55])
    c = np.array([0,1/4,3/8,12/13,1,1/2])
    d = np.array([1/360,0,-128/4275,-2197/75240,1/50,2/55])
    k = np.zeros_like(c)
    for i in range(6):
        if i == 0:
            k[i] = delta_t * func(y0,t0)
        else:
            k[i] = delta_t * func(y0 + (np.sum(a[i,:] * k[:])),t0 + c[i] * delta_t)
    TE = np.abs(np.sum(k * d))
    print(TE)
    if tolerance <= TE:
        delta_t = 0.9 * delta_t * (tolerance / TE) ** (1/5)
        tf = t0 + delta_t
        for i in range(6):
            if i == 0:
                k[i] = delta_t * func(y0,t0)
            else:
                k[i] = delta_t * func(y0 + (np.sum(a[i,:] * k[:])),t0 + c[i] * delta_t)
    yf = y0 + np.sum(b5 * k)
    return yf, tf
