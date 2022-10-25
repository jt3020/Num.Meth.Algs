import numpy as np
import sys

# 2D Transfinite Interpolation
# Performed on eight boundaries
# Bottom, Top, Left, Right
def Transfinite_Interpolation(xb,yb,xt,yt,xl,yl,xr,yr,x,y):
    
    for j in range(2,n):
        s = (j-1)*dy

        for i in range(2,m):
            r = (i-1)*dx

            x(i,j) = (1.-s)*xb(i) + s*xt(i) + (1.-r)*xl(j)   \
                           + r*xr(j) - ( r*s*xr(n+1) +       \
                           (1.-r)*s*xl(n+1) + r*(1.-s)*xr(1) \
                           + (1.-r)*(1.-s)*xl(1) )

            y(i,j) = (1.-s)*yb(i) + s*yt(i) + (1.-r)*yl(j)   \
                           + r*yr(j) - ( r*s*yr(n+1) +       \
                           (1.-r)*s*yl(n+1) + r*(1.-s)*yr(1) \
                           + (1.-r)*(1.-s)*yl(1) )

