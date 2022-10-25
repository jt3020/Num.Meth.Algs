import numpy as np

def Backward_Euler(func,y0,t0,tf,tolerance,delta):
    """Calculate the solution of func(y,t) at tf with initial condition y0 at t0"""
    delta_t = tf - t0
    print("delta t:", delta_t)
    yf = y0
    print("yf:", yf)
    print("func:", func(yf,tf))
    residual = yf - y0 - delta_t * func(yf,tf)
    print("res :",residual)
    while np.absolute(residual) > tolerance:
        
        m = ((yf + delta - y0 - delta_t * func(yf + delta,tf)) - (yf - y0 - delta_t * func(yf,tf))) / delta
        yf = yf - (yf - y0 - delta_t * func(yf,tf)) / m
        residual = yf - y0 - delta_t * func(yf,tf)
        print(residual)
    return yf
