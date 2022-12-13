import numpy as np
import sys
import matplotlib.pyplot as plt

# 2D Transfinite Interpolation
# Performed on eight boundaries
# Bottom, Top, Left, Right)

def Transfinite_Interpolation(xb,yb,xt,yt,xl,yl,xr,yr,x,y):

    n = len(xb)-1
    m = len(xl)-1
    dx = 1/n 
    dy = 1/m
    
    for j in range(1,m):
        Eta = (j)*dy
        for i in range(1,n):
            Xi = (i)*dx
            x[j,i] = (1.-Eta)*xb[i] + Eta*xt[i] + (1.-Xi)*xl[j]  + Xi*xr[j] \
                    - ( Xi*Eta*xr[m] + (1.-Xi)*Eta*xl[m] + Xi*(1.-Eta)*xr[0] + (1.-Xi)*(1.-Eta)*xl[0] )

            y[j,i] = (1.-Eta)*yb[i] + Eta*yt[i] + (1.-Xi)*yl[j]  + Xi*yr[j] \
                    - ( Xi*Eta*yr[m] + (1.-Xi)*Eta*yl[m] + Xi*(1.-Eta)*yr[0] + (1.-Xi)*(1.-Eta)*yl[0] )

n = 20
m = 10
dx = 1/(n-1)
dy = 1/(m-1)

x = np.zeros((m,n))
y = np.zeros((m,n))

xb = np.zeros(n)
yb = np.zeros(n)

xt = np.zeros(n)
yt = np.zeros(n)

xl = np.zeros(m)
yl = np.zeros(m)

xr = np.zeros(m)
yr = np.zeros(m)

for i in range(n):
    xb[i] = i/(n-1)
    yb[i] = 0.0 + 0.5*np.sin(np.pi*xb[i])
    xt[i] = i/(n-1)
    yt[i] = 1.0 + 0.5*np.sin(np.pi*xb[i])

for i in range(m):
    xl[i] = 0.0
    yl[i] = i/(m-1)
    xr[i] = 1.0
    yr[i] = i/(m-1)

x[0,:] = xb
x[m-1,:] = xt
x[:,0] = xl 
x[:,n-1] = xr

y[0,:] = yb
y[m-1,:] = yt
y[:,0] = yl 
y[:,n-1] = yr

Transfinite_Interpolation(xb,yb,xt,yt,xl,yl,xr,yr,x,y)

for i in range(n):
    plt.plot(x[:,i],y[:,i],color='black')
for i in range(m):
    plt.plot(x[i,:],y[i,:],color='black')
plt.show()