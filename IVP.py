import numpy as np
import matplotlib.pyplot as plt
import ODE_solvers as solver_single
import ODEs_solvers as solver_system

def single_solve(func,y0,t0,tf,timestep,solver):
    """Soves an ODE with initial conditions t0 at t0"""
    y = np.array([y0])
    t = np.array([t0])
    while t[-1] < tf:
        if solver == "Backward Euler":
            y = np.append(y,[solver_single.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "Crank Nicolson":
            y = np.append(y,[solver_single.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "RK45":
            y = np.append(y,[solver_single.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[1]],axis=0)
        elif solver == "A_B_M":
            y = np.append(y,[solver_single.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[0]],axis=0)
            t = np.append(t,[solver_single.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[1]],axis=0)
        else:
            print("select a solver")
            break
    return y,t

def system_solve(func,y0,t0,tf,timestep,solver):
    """Soves an ODE with initial conditions t0 at t0"""
    y = np.array([y0])
    t = np.array([t0])
    while t[-1] < tf:
        if solver == "Backward Euler":
            y = np.append(y,[solver_single.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.Backward_Euler(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "Crank Nicolson":
            y = np.append(y,[solver_single.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.Crank_Nicolson(func,y[-1],t[-1],t[-1]+timestep,0.001,0.001)[1]],axis=0)
        elif solver == "RK45":
            y = np.append(y,[solver_single.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[0]],axis=0)
            t = np.append(t,[solver_single.RK45(func,y[-1],t[-1],t[-1]+timestep,0.001)[1]],axis=0)
        elif solver == "A_B_M":
            y = np.append(y,[solver_single.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[0]],axis=0)
            t = np.append(t,[solver_single.A_B_M(func,y[-1],t[-1],t[-1]+timestep)[1]],axis=0)
        else:
            print("select a solver")
            break
    return y,t

def test(y,t):
    res = np.zeros(2)
    res[0] = y[0] - y[0] * y[1]
    res[1] = -y[1] + y[0] * y[1]
    return res

print(solver_system.A_B_M(test,np.array([1]),0,0.1))

y,t = single_solve(test,0,0,10,2,"A_B_M")

plt.figure()
plt.plot(t,y)
plt.plot(t,np.sin(t))
plt.show()
