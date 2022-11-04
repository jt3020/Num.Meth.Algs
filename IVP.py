import numpy as np
import matplotlib.pyplot as plt
import ODE_solvers as solvers

def solve(func,y0,t0,tf,timestep,solver):
    """Soves an ODE with initial conditions t0 at t0"""
    y = np.array([y0])
    t = np.array([t0])
    while t[-1] < tf:
        if solver == "Backward Euler":
            y = np.append(y,[solvers.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solvers.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "Crank Nicolson":
            y = np.append(y,[solvers.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solvers.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "RK45":
            y = np.append(y,[solvers.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[0]],axis=0)
            t = np.append(t,[solvers.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[1]],axis=0)
        elif solver == "A_B_M":
            y = np.append(y,[solvers.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[0]],axis=0)
            t = np.append(t,[solvers.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[1]],axis=0)
        else:
            print("select a solver")
            break
    return y,t

def test(y,t):
    return -y

y,t=solve(test,1,0,1,0.1,'RK45')

